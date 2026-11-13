// rl_analysis <output-dir>
//
// Runs every analysis behind docs/RESULTS.md and writes one CSV per table into <output-dir>.
// tools/make_results.py turns those into the tables and plots. Everything is deterministic (fixed
// seeds), so the committed numbers can be regenerated exactly; only the timings depend on the machine.
//
// Units: prices in basis points of notional (1 bp = 0.01% of notional), vols in bp of normal vol.
#include "rl/bootstrap/global_bootstrap.hpp"
#include "rl/bootstrap/jacobian.hpp"
#include "rl/bootstrap/report.hpp"
#include "rl/core/error.hpp"
#include "rl/instruments/builder.hpp"
#include "rl/linear/forward.hpp"
#include "rl/linear/leg_pv.hpp"
#include "rl/models/g2pp_calibration.hpp"
#include "rl/models/g2pp_swaption.hpp"
#include "rl/models/lgm_swaption.hpp"
#include "rl/math/rng.hpp"
#include "rl/models/usd_market.hpp"
#include "rl/pricers/bermudan_g2pp.hpp"
#include "rl/pricers/bermudan_mc.hpp"
#include "rl/pricers/bermudan_pde.hpp"
#include "rl/pricers/cms_spread.hpp"
#include "rl/pricers/exposure.hpp"
#include "rl/pricers/risk.hpp"
#include "rl/vol/arbitrage.hpp"
#include "rl/xva/bermudan_exposure.hpp"
#include "rl/xva/credit.hpp"
#include "rl/xva/cva_risk.hpp"
#include "rl/xva/exposure_engine.hpp"
#include "rl/xva/wrong_way.hpp"
#include "rl/xva/xva.hpp"
#include "rl/vol/sabr.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <functional>
#include <iostream>
#include <string>
#include <tuple>
#include <vector>

using namespace rl;
using market::CurveInterpolation;
using pricers::BermudanProblem;
using pricers::CalibratedBermudan;

namespace {

const double kBp = 1e4;
const double kKappa = 0.03; // reference mean reversion

struct Csv {
    std::FILE* f;
    Csv(const std::filesystem::path& dir, const char* name, const char* header) : f(std::fopen((dir / name).string().c_str(), "w")) {
        if (!f) throw InputError(std::string("cannot write ") + name);
        std::fprintf(f, "%s\n", header);
    }
    ~Csv() { std::fclose(f); }
    Csv(const Csv&) = delete;
    Csv& operator=(const Csv&) = delete;
};

/// Milliseconds per call, repeating until at least 0.2 s has elapsed.
double timeMs(const std::function<void()>& fn) {
    using clock = std::chrono::steady_clock;
    int reps = 0;
    const auto t0 = clock::now();
    double elapsed = 0.0;
    do {
        fn();
        ++reps;
        elapsed = std::chrono::duration<double, std::milli>(clock::now() - t0).count();
    } while (elapsed < 200.0);
    return elapsed / reps;
}

instruments::Swap swap(const models::UsdMarket& m, int years, double rate, bool payer = true, double notional = 1.0) {
    return instruments::OisSwapBuilder::usdSofr(m.asOf).id("swap").tenor(core::Tenor::years(years)).fixedRate(rate).notional(notional).payFixed(payer).build();
}

pricers::PdeOptions pdeGrid(int n, bool rannacher = true) {
    pricers::PdeOptions o;
    o.spaceSteps = n;
    o.timeSteps = n;
    o.rannacher = rannacher;
    return o;
}

// ---------------------------------------------------------------------------------------------------
// Curves
// ---------------------------------------------------------------------------------------------------

void curves(const std::filesystem::path& out) {
    {
        Csv solver(out, "curve_solver.csv", "interpolation,method,iterations,max_residual_bp,worst_instrument,time_ms");
        Csv rep(out, "curve_reprice.csv", "interpolation,id,type,curve,quote,model,residual_bp");
        for (auto interp : {CurveInterpolation::LinearLogDf, CurveInterpolation::MonotoneConvex}) {
            const models::UsdMarket m = models::loadUsdMarket(interp);
            const std::string name(market::toString(interp));
            auto worst = [&](const bootstrap::BootstrapResult& r) {
                Eigen::Index i = 0;
                r.scaledResiduals.cwiseAbs().maxCoeff(&i);
                return m.problem.instrumentIds()[static_cast<std::size_t>(i)];
            };
            const bootstrap::BootstrapResult glob = bootstrap::solve(m.problem), seq = bootstrap::solveSequential(m.problem);
            bootstrap::BootstrapOptions fd;
            fd.analyticJacobian = false;
            const bootstrap::BootstrapResult globFd = bootstrap::solve(m.problem, fd);
            std::fprintf(solver.f, "%s,global Newton (adjoint Jacobian),%d,%.3e,%s,%.3f\n", name.c_str(), glob.iterations, glob.maxAbsBp, worst(glob).c_str(),
                         timeMs([&] { bootstrap::solve(m.problem); }));
            std::fprintf(solver.f, "%s,global Newton (finite-difference Jacobian),%d,%.3e,%s,%.3f\n", name.c_str(), globFd.iterations, globFd.maxAbsBp,
                         worst(globFd).c_str(), timeMs([&] { bootstrap::solve(m.problem, fd); }));
            std::fprintf(solver.f, "%s,sequential (one pass curve by curve),%d,%.3e,%s,%.3f\n", name.c_str(), seq.iterations, seq.maxAbsBp, worst(seq).c_str(),
                         timeMs([&] { bootstrap::solveSequential(m.problem); }));
            for (const auto& row : bootstrap::reprice(m.problem, m.z).rows)
                std::fprintf(rep.f, "%s,%s,%s,%s,%.6f,%.6f,%.3e\n", name.c_str(), row.id.c_str(), std::string(market::toString(row.type)).c_str(), row.curve.c_str(),
                             row.quote, row.model, row.residualBp);
        }
    }
    {
        // Instantaneous-style forwards (1-week continuous) under the two interpolations.
        Csv fwd(out, "curve_forwards.csv", "t,linear_log_df,monotone_convex");
        const models::UsdMarket lin = models::loadUsdMarket(CurveInterpolation::LinearLogDf), mc = models::loadUsdMarket(CurveInterpolation::MonotoneConvex);
        const double w = 7.0 / 365.0;
        for (double t = 0.02; t < 30.0; t += 0.02)
            std::fprintf(fwd.f, "%.4f,%.6f,%.6f\n", t, 100.0 * linear::continuousForward(lin.curves.discount(), t, t + w),
                         100.0 * linear::continuousForward(mc.curves.discount(), t, t + w));
    }
}

// ---------------------------------------------------------------------------------------------------
// Swap risk: adjoint against bump
// ---------------------------------------------------------------------------------------------------

void swapRisk(const std::filesystem::path& out, const models::UsdMarket& m) {
    // A book of 200 swaps: maturities 1..30y, alternating direction, strikes around the market.
    std::vector<instruments::Swap> book;
    for (int i = 0; i < 200; ++i) book.push_back(swap(m, 1 + (i * 7) % 30, 0.025 + 0.0001 * (i % 150), i % 3 != 0, 1e6 * (1 + i % 5)));
    const auto aad = pricers::swapParDeltaAdjoint(m, book), bump = pricers::swapParDeltaBump(m, book), small = pricers::swapParDeltaBump(m, book, 0.01);
    Csv rows(out, "swap_risk.csv", "quote,adjoint,bump_1bp,bump_0.01bp");
    double largest = 0.0, worst = 0.0;
    for (std::size_t i = 0; i < aad.size(); ++i) {
        std::fprintf(rows.f, "%s,%.4f,%.4f,%.4f\n", aad[i].id.c_str(), aad[i].value, bump[i].value, small[i].value);
        largest = std::max(largest, std::fabs(aad[i].value));
        worst = std::max(worst, std::fabs(aad[i].value - small[i].value));
    }
    Csv t(out, "swap_risk_timing.csv", "method,swaps,quotes,time_ms,max_difference_rel");
    const double tPv = timeMs([&] {
        double v = 0.0;
        for (const auto& s : book) v += linear::pv(s, m.curves, m.fixings);
        if (v == 1e300) std::puts("");
    });
    const double tAad = timeMs([&] { pricers::swapParDeltaAdjoint(m, book); });
    const double tBump = timeMs([&] { pricers::swapParDeltaBump(m, book); });
    std::fprintf(t.f, "one valuation of the book,%zu,%zu,%.3f,\n", book.size(), aad.size(), tPv);
    std::fprintf(t.f, "adjoint + implicit-function theorem,%zu,%zu,%.3f,%.2e\n", book.size(), aad.size(), tAad, worst / largest);
    std::fprintf(t.f, "bump and re-bootstrap each quote,%zu,%zu,%.3f,\n", book.size(), aad.size(), tBump);
}

// ---------------------------------------------------------------------------------------------------
// Vol cube
// ---------------------------------------------------------------------------------------------------

void volCube(const std::filesystem::path& out, const models::UsdMarket& m) {
    Csv fit(out, "sabr_fit.csv", "node,expiry,tenor,forward_pct,atm_bp,alpha,rho,nu,max_residual_bp,iterations,ok,first_negative_density_strike_pct");
    const auto& cube = m.cube;
    for (std::size_t i = 0; i < cube.expiryCount(); ++i)
        for (std::size_t j = 0; j < cube.tenorCount(); ++j) {
            const std::size_t k = i * cube.tenorCount() + j;
            const auto& node = cube.node(i, j);
            const auto& f = m.cubeFit.fits[k];
            const double e = cube.expiries()[i];
            const auto range = vol::diagnosticStrikeRange(f.params, node.forward, e, node.atmNormalVol);
            const vol::ArbitrageCheck chk = vol::checkSabrSmile(f.params, node.forward, e, range.first, range.second, 400);
            std::fprintf(fit.f, "%s,%.4f,%.1f,%.4f,%.2f,%.6f,%.4f,%.4f,%.4f,%d,%d,%s\n", m.cubeFit.nodeIds[k].c_str(), e, cube.tenors()[j], node.forward * 100.0,
                         node.atmNormalVol * kBp, f.params.alpha, f.params.rho, f.params.nu, f.maxAbsBp, f.iterations, f.ok ? 1 : 0,
                         chk.firstNegativeDensityBelowForward ? std::to_string(*chk.firstNegativeDensityBelowForward * 100.0).c_str() : "");
        }
    // Two smiles with their densities: a benign node and the stressed 10Y x 30Y node.
    Csv smile(out, "sabr_smile.csv", "node,strike_pct,normal_vol_bp,density");
    for (const char* id : {"5Yx10Y", "10Yx30Y"}) {
        for (std::size_t k = 0; k < m.cubeFit.nodeIds.size(); ++k) {
            if (m.cubeFit.nodeIds[k] != id) continue;
            const auto& node = cube.nodes()[k];
            const double e = cube.expiries()[k / cube.tenorCount()];
            const auto strikes = vol::strikeGrid(std::max(node.forward - 0.03, 0.002), node.forward + 0.03, 241);
            for (const auto& p : vol::riskNeutralDensity(m.cubeFit.fits[k].params, node.forward, e, strikes))
                std::fprintf(smile.f, "%s,%.4f,%.4f,%.6f\n", id, p.strike * 100.0, p.normalVol * kBp, p.density);
        }
    }
}

// ---------------------------------------------------------------------------------------------------
// LGM calibration
// ---------------------------------------------------------------------------------------------------

void lgm(const std::filesystem::path& out, const models::UsdMarket& m) {
    Csv cal(out, "lgm_calibration.csv", "kappa,expiry,label,tenor,forward_pct,strike_pct,market_vol_bp,model_vol_bp,error_bp,zeta,hw_sigma_bp");
    for (double kappa : {0.0, 0.03, 0.10}) {
        const CalibratedBermudan b = pricers::calibrateBermudan(m, swap(m, 10, 0.035), 1.0, 0, kappa);
        for (const auto& r : b.calibration.rows)
            std::fprintf(cal.f, "%.2f,%.4f,%s,%.2f,%.4f,%.4f,%.4f,%.4f,%.2e,%.6e,%.3f\n", kappa, r.expiry, r.label.c_str(), r.tenor, r.forward * 100.0,
                         r.strike * 100.0, r.marketVolBp, r.modelVolBp, r.errorBp, r.zeta, r.hullWhiteSigmaBp);
    }
}

// ---------------------------------------------------------------------------------------------------
// Bermudan
// ---------------------------------------------------------------------------------------------------

void bermudan(const std::filesystem::path& out, const models::UsdMarket& m) {
    {
        Csv t(out, "bermudan_three_ways.csv",
              "trade,max_european_bp,pde_bp,lsmc_lower_bp,lsmc_se_bp,upper_bp,upper_se_bp,gap_bp,gap_se_bp,pde_ms,lsmc_ms,upper_ms,exercise_probability");
        struct Trade {
            const char* name;
            double strike;
            bool payer;
        };
        for (const Trade& tr : {Trade{"10nc1 payer 2.50% (ITM)", 0.025, true}, Trade{"10nc1 payer 3.50% (near ATM)", 0.035, true}, Trade{"10nc1 payer 4.50% (OTM)", 0.045, true},
                                Trade{"10nc1 receiver 3.50%", 0.035, false}}) {
            const CalibratedBermudan b = pricers::calibrateBermudan(m, swap(m, 10, tr.strike, tr.payer), 1.0, 0, kKappa);
            const auto european = b.problem.europeanValues();
            // Independent random numbers per trade, so the four rows are four independent experiments.
            pricers::LsmcOptions lo;
            lo.pricingPaths = 4000000;
            lo.seed = math::SeedPolicy::derive(20260904, static_cast<std::uint64_t>(tr.strike * 1e4), tr.payer ? 1 : 2);
            pricers::UpperBoundOptions uo;
            uo.seed = lo.seed + 1;
            const double pde = pricers::priceBermudanPde(b.problem, pdeGrid(1600)).value;
            const pricers::LsmcResult low = pricers::priceBermudanLsmc(b.problem, lo);
            const pricers::UpperBoundResult up = pricers::priceBermudanUpperBound(b.problem, low, uo);
            double prob = 0.0;
            for (double q : low.exerciseProbability) prob += q;
            std::fprintf(t.f, "%s,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.1f,%.1f,%.1f,%.3f\n", tr.name, *std::max_element(european.begin(), european.end()) * kBp,
                         pde * kBp, low.lower * kBp, low.standardError * kBp, up.upper * kBp, up.standardError * kBp, up.gap * kBp, up.gapStandardError * kBp,
                         timeMs([&] { pricers::priceBermudanPde(b.problem, pdeGrid(400)); }), timeMs([&] { pricers::priceBermudanLsmc(b.problem, lo); }),
                         timeMs([&] { pricers::priceBermudanUpperBound(b.problem, low, uo); }), prob);
        }
    }
    const CalibratedBermudan ref = pricers::calibrateBermudan(m, swap(m, 10, 0.035), 1.0, 0, kKappa);
    const double exact = pricers::priceBermudanPde(ref.problem, pdeGrid(6400, false)).value;
    {
        Csv c(out, "pde_convergence.csv", "steps,error_rannacher_bp,error_cn_bp,time_ms");
        for (int n : {50, 100, 200, 400, 800, 1600, 3200})
            std::fprintf(c.f, "%d,%.3e,%.3e,%.3f\n", n, std::fabs(pricers::priceBermudanPde(ref.problem, pdeGrid(n, true)).value - exact) * kBp,
                         std::fabs(pricers::priceBermudanPde(ref.problem, pdeGrid(n, false)).value - exact) * kBp,
                         timeMs([&] { pricers::priceBermudanPde(ref.problem, pdeGrid(n, true)); }));
        Csv ct(out, "pde_time_convergence.csv", "time_steps,error_rannacher_bp,error_cn_bp");
        for (int n : {25, 50, 100, 200, 400, 800}) {
            pricers::PdeOptions a = pdeGrid(3200, true), b = pdeGrid(3200, false);
            a.timeSteps = b.timeSteps = n;
            std::fprintf(ct.f, "%d,%.3e,%.3e\n", n, std::fabs(pricers::priceBermudanPde(ref.problem, a).value - exact) * kBp,
                         std::fabs(pricers::priceBermudanPde(ref.problem, b).value - exact) * kBp);
        }
        // Gamma in the state today, on a fine space grid with few time steps.
        Csv g(out, "pde_gamma.csv", "x,gamma_rannacher,gamma_cn");
        pricers::PdeOptions a = pdeGrid(3200, true), b = pdeGrid(3200, false);
        a.timeSteps = b.timeSteps = 36;
        const pricers::PdeResult ra = pricers::priceBermudanPde(ref.problem, a), rb = pricers::priceBermudanPde(ref.problem, b);
        const double h = ra.grid[1] - ra.grid[0];
        const std::size_t c0 = ra.grid.size() / 2;
        for (std::size_t j = c0 - 500; j <= c0 + 500; j += 2)
            std::fprintf(g.f, "%.6f,%.6f,%.6f\n", ra.grid[j], (ra.values[j - 1] - 2 * ra.values[j] + ra.values[j + 1]) / (h * h),
                         (rb.values[j - 1] - 2 * rb.values[j] + rb.values[j + 1]) / (h * h));
    }
    {
        // Regression basis: quality of the exercise rule, measured by both bounds.
        Csv c(out, "lsmc_basis.csv", "basis,degree,lower_bp,lower_se_bp,gap_bp,gap_se_bp,pde_bp");
        for (auto basis : {pricers::RegressionBasis::State, pricers::RegressionBasis::ExerciseValue})
            for (int degree : {0, 1, 2, 3, 4, 6}) {
                pricers::LsmcOptions lo;
                lo.basis = basis;
                lo.degree = degree;
                lo.pricingPaths = 400000;
                const auto low = pricers::priceBermudanLsmc(ref.problem, lo);
                const auto up = pricers::priceBermudanUpperBound(ref.problem, low);
                std::fprintf(c.f, "%s,%d,%.3f,%.3f,%.3f,%.3f,%.3f\n", basis == pricers::RegressionBasis::State ? "state" : "exercise value", degree, low.lower * kBp,
                             low.standardError * kBp, up.gap * kBp, up.gapStandardError * kBp, exact * kBp);
            }
    }
    {
        // In-sample bias against the number of training paths (degree 6, averaged over 200 seeds).
        Csv c(out, "lsmc_bias.csv", "training_paths,in_sample_bp,out_of_sample_bp,pde_bp");
        for (std::size_t n : {250u, 500u, 1000u, 2000u, 4000u, 8000u, 16000u}) {
            double in = 0.0, outv = 0.0;
            const int trials = 200;
            for (int k = 0; k < trials; ++k) {
                pricers::LsmcOptions lo;
                lo.trainingPaths = n;
                lo.pricingPaths = 10000;
                lo.degree = 6;
                lo.seed = 500 + static_cast<std::uint64_t>(k);
                const auto r = pricers::priceBermudanLsmc(ref.problem, lo);
                in += r.inSample / trials;
                outv += r.lower / trials;
            }
            std::fprintf(c.f, "%zu,%.3f,%.3f,%.3f\n", n, in * kBp, outv * kBp, exact * kBp);
        }
    }
    {
        // Monte Carlo error against the number of pricing paths.
        Csv c(out, "lsmc_paths.csv", "pricing_paths,standard_error_bp,standard_error_no_antithetic_bp,time_ms");
        for (std::size_t n : {1000u, 4000u, 16000u, 64000u, 256000u, 1024000u}) {
            pricers::LsmcOptions lo, plain;
            lo.pricingPaths = plain.pricingPaths = n;
            plain.antithetic = false;
            std::fprintf(c.f, "%zu,%.4f,%.4f,%.2f\n", n, pricers::priceBermudanLsmc(ref.problem, lo).standardError * kBp,
                         pricers::priceBermudanLsmc(ref.problem, plain).standardError * kBp, timeMs([&] { pricers::priceBermudanLsmc(ref.problem, lo); }));
        }
    }
    {
        // Mean reversion.
        Csv c(out, "bermudan_kappa.csv", "kappa,bermudan_bp,max_european_bp,switch_value_bp,hw_sigma_first_bp,hw_sigma_last_bp");
        for (double kappa : {-0.02, -0.01, 0.0, 0.01, 0.02, 0.03, 0.05, 0.075, 0.10, 0.15}) {
            const CalibratedBermudan b = pricers::calibrateBermudan(m, swap(m, 10, 0.035), 1.0, 0, kappa);
            const auto european = b.problem.europeanValues();
            const double best = *std::max_element(european.begin(), european.end()), v = pricers::priceBermudanPde(b.problem, pdeGrid(800)).value;
            std::fprintf(c.f, "%.3f,%.3f,%.3f,%.3f,%.2f,%.2f\n", kappa, v * kBp, best * kBp, (v - best) * kBp, b.calibration.rows.front().hullWhiteSigmaBp,
                         b.calibration.rows.back().hullWhiteSigmaBp);
        }
    }
    {
        // Exercise rights: value against how many of the annual dates are kept (last dates first).
        Csv c(out, "bermudan_dates.csv", "exercise_dates,first_date,value_bp");
        std::vector<models::ExerciseRight> keep;
        for (std::size_t i = ref.problem.size(); i-- > 0;) {
            keep.insert(keep.begin(), ref.problem.rights()[i]);
            std::fprintf(c.f, "%zu,%s,%.3f\n", keep.size(), keep.front().label.c_str(), pricers::priceBermudanPde(BermudanProblem(ref.problem.model(), keep), pdeGrid(800)).value * kBp);
        }
        // Exercise boundary and probabilities.
        Csv e(out, "bermudan_boundary.csv", "date,expiry,forward_pct,european_bp,boundary_rate_pct,exercise_probability");
        const pricers::PdeResult pde = pricers::priceBermudanPde(ref.problem, pdeGrid(1600));
        pricers::LsmcOptions lo;
        lo.pricingPaths = 400000;
        const pricers::LsmcResult low = pricers::priceBermudanLsmc(ref.problem, lo);
        const auto european = ref.problem.europeanValues();
        for (std::size_t i = 0; i < ref.problem.size(); ++i)
            std::fprintf(e.f, "%s,%.3f,%.4f,%.3f,%.4f,%.4f\n", ref.problem.rights()[i].label.c_str(), ref.problem.rights()[i].time,
                         ref.problem.rights()[i].forwardRate() * 100.0, european[i] * kBp, pde.boundaryRate[i] * 100.0, low.exerciseProbability[i]);
    }
    {
        // Risk, notional 1,000,000.
        const instruments::Swap s = swap(m, 10, 0.035, true, 1e6);
        auto price = [&](const models::UsdMarket& mk) { return pricers::priceBermudanPde(pricers::calibrateBermudan(mk, s, 1.0, 0, kKappa).problem, pdeGrid(400)).value; };
        const auto bermudanDelta = pricers::parDeltaBump(m, price);
        // For comparison: the swap from the first call date (what the option becomes if exercised at once).
        const auto firstCall = pricers::calibrateBermudan(m, s, 1.0, 0, kKappa).problem.rights().front().date;
        const instruments::Swap fwd = instruments::OisSwapBuilder::usdSofr(m.asOf).id("fwd").effective(firstCall).maturity(m.asOf.addYears(10).addDays(4)).fixedRate(0.035).notional(1e6).build();
        const auto swapDelta = pricers::swapParDeltaAdjoint(m, {fwd});
        Csv d(out, "bermudan_delta.csv", "quote,bermudan,forward_swap");
        for (std::size_t i = 0; i < bermudanDelta.size(); ++i)
            if (std::fabs(bermudanDelta[i].value) > 0.005 || std::fabs(swapDelta[i].value) > 0.005)
                std::fprintf(d.f, "%s,%.3f,%.3f\n", bermudanDelta[i].id.c_str(), bermudanDelta[i].value, swapDelta[i].value);
        Csv v(out, "bermudan_vega.csv", "expiry_date,bermudan_vega,european_vega");
        const auto vega = pricers::bermudanVegaByExpiry(m, s, 1.0, 0, kKappa, [](const BermudanProblem& p) { return pricers::priceBermudanPde(p, pdeGrid(400)).value; });
        const CalibratedBermudan b = pricers::calibrateBermudan(m, s, 1.0, 0, kKappa);
        for (std::size_t i = 0; i < vega.size(); ++i) {
            const auto& r = b.problem.rights()[i];
            const double vol0 = b.calibration.rows[i].marketVolBp / kBp;
            const double europeanVega = models::marketSwaptionPrice(r, vol0 + 1e-4) - models::marketSwaptionPrice(r, vol0);
            std::fprintf(v.f, "%s,%.3f,%.3f\n", vega[i].id.c_str(), vega[i].value, europeanVega);
        }
    }
}

// ---------------------------------------------------------------------------------------------------
// Two factors (G2++)
// ---------------------------------------------------------------------------------------------------

void twoFactor(const std::filesystem::path& out, const models::UsdMarket& m) {
    using models::G2pp;
    const auto matrix = models::makeAtmMatrix(m.asOf, m.curves, m.cube, {1, 2, 3, 5, 7, 10}, {1, 2, 5, 10, 20});

    // --- Shape from the ATM matrix: one factor, then two factors from several starting points -------
    Csv par(out, "g2pp_params.csv", "model,start,kappa1,kappa2,sigma_ratio,rho,rms_bp,max_bp,iterations,status,active_bounds,corr_1y_10y,corr_2y_30y,fit_ms");
    auto report = [&](const char* name, const char* start, const models::G2ppFit& f, double ms) {
        std::string bounds;
        for (const std::string& b : f.activeBounds) bounds += (bounds.empty() ? "" : " ") + b;
        std::fprintf(par.f, "%s,%s,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%d,%s,%s,%.4f,%.4f,%.0f\n", name, start, f.model.kappa1(), f.model.kappa2(),
                     f.model.sigma2() / f.model.sigma1(), f.model.rho(), f.rmsErrorBp, f.maxAbsErrorBp, f.iterations, math::toString(f.status), bounds.c_str(),
                     f.model.forwardRateCorrelation(1.0, 10.0), f.model.forwardRateCorrelation(2.0, 30.0), ms);
    };
    auto timedFit = [&](const models::G2ppFitOptions& o, double& ms) {
        const auto t0 = std::chrono::steady_clock::now();
        models::G2ppFit f = models::fitG2pp(matrix, o);
        ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
        return f;
    };
    double ms = 0.0;
    models::G2ppFitOptions one;
    one.factors = 1;
    one.kappa1 = 0.05;
    const models::G2ppFit fit1 = timedFit(one, ms);
    report("1F", "kappa 0.05", fit1, ms);
    std::vector<models::G2ppFit> fits2;
    const struct {
        double kappa1, rho;
        const char* label;
    } starts[] = {{0.1, -0.9, "kappa1 0.1 rho -0.9"}, {0.1, 0.5, "kappa1 0.1 rho 0.5"}, {0.6, -0.9, "kappa1 0.6 rho -0.9"}, {0.6, 0.5, "kappa1 0.6 rho 0.5"}};
    for (const auto& s : starts) {
        models::G2ppFitOptions o;
        o.kappa1 = s.kappa1;
        o.rho = s.rho;
        fits2.push_back(timedFit(o, ms));
        report("2F", s.label, fits2.back(), ms);
    }
    const models::G2ppFit& fit2 = *std::min_element(fits2.begin(), fits2.end(), [](const auto& a, const auto& b) { return a.rmsErrorBp < b.rmsErrorBp; });
    {
        Csv c(out, "g2pp_fit.csv", "label,expiry,tenor,market_vol_bp,error_1f_bp,error_2f_bp");
        for (std::size_t i = 0; i < matrix.size(); ++i)
            std::fprintf(c.f, "%s,%.3f,%.3f,%.3f,%.3f,%.3f\n", fit1.rows[i].label.c_str(), fit1.rows[i].expiry, fit1.rows[i].tenor, fit1.rows[i].marketVolBp,
                         fit1.rows[i].errorBp, fit2.rows[i].errorBp);
    }

    // --- The swaption formula: accuracy and cost -----------------------------------------------------
    {
        Csv c(out, "g2pp_swaption.csv", "quadrature_points,price_bp,error_vs_2d_integration_bp,us_per_price");
        const G2pp model(0.40, 0.03, 0.012, 0.009, -0.6);
        const models::ExerciseRight& r = matrix[18].right; // 5Y x 10Y
        const double ref = models::g2ppOptionValueByIntegration(model, r.time, r.flows, 2401);
        for (int n : {2, 4, 8, 16, 24}) {
            const double v = models::g2ppOptionValue(model, r.time, r.flows, n);
            const double us = 1000.0 * timeMs([&] { (void)models::g2ppOptionValue(model, r.time, r.flows, n); });
            std::fprintf(c.f, "%d,%.6f,%.2e,%.1f\n", n, v * kBp, (v - ref) * kBp, us);
        }
    }

    // --- The Bermudan: same Europeans, different correlation ---------------------------------------
    // Every row reprices the nine co-terminal swaptions of the 10nc1 payer at 3.50%. What changes is
    // the shape: with the second factor off the model is the one-factor model; with it on, rates on
    // different exercise dates decorrelate by an amount set by rho.
    Csv b(out, "g2pp_bermudan.csv", "case,kappa1,kappa2,sigma_ratio,rho,corr_1y_10y,strip_error_bp,price_bp,se_bp,in_sample_bp,best_european_bp,method");
    const instruments::Swap trade = swap(m, 10, 0.035);
    const CalibratedBermudan ref = pricers::calibrateBermudan(m, trade, 1.0, 0, kKappa);
    const std::vector<double> europeans = ref.problem.europeanValues();
    const double best = *std::max_element(europeans.begin(), europeans.end());
    std::fprintf(b.f, "one factor,%.4f,,,,1.0000,%.1e,%.3f,,,%.3f,PDE\n", kKappa, ref.calibration.maxAbsErrorBp, pricers::priceBermudanPde(ref.problem, pdeGrid(800)).value * kBp,
                 best * kBp);
    pricers::G2ppLsmcOptions mc;
    mc.trainingPaths = 50000;
    mc.pricingPaths = 1000000;
    auto row = [&](const char* name, const G2pp& shape) {
        const pricers::CalibratedBermudanG2pp cal = pricers::calibrateBermudanG2pp(m, trade, 1.0, 0, shape);
        const pricers::G2ppLsmcResult r = pricers::priceBermudanLsmcG2pp(cal.calibration.model, cal.rights, mc);
        std::fprintf(b.f, "%s,%.4f,%.4f,%.4f,%.4f,%.4f,%.1e,%.3f,%.3f,%.3f,%.3f,LSMC\n", name, shape.kappa1(), shape.kappa2(), shape.sigma2() / shape.sigma1(), shape.rho(),
                     shape.forwardRateCorrelation(1.0, 10.0), cal.calibration.maxAbsErrorBp, r.lower * kBp, r.standardError * kBp, r.inSample * kBp, best * kBp);
    };
    row("second factor off", G2pp(kKappa, 1.0, 0.01, 0.0, 0.0));
    for (double rho : {0.9, 0.5, 0.0, -0.5, -0.9}) row("fast second factor", G2pp(kKappa, 0.40, 0.01, 0.01, rho));
    row("fitted to the ATM matrix", G2pp(fit2.model.kappa1(), fit2.model.kappa2(), fit2.model.sigma1(), fit2.model.sigma2(), fit2.model.rho()));

    // --- A CMS spread option: the product that prices the correlation ------------------------------
    // 5-year expiry on the 10-year minus the 2-year swap rate. Every model below reprices the 5Yx2Y
    // and 5Yx10Y swaptions exactly, so both rates' volatilities agree; what differs is how the two
    // rates move together.
    {
        const auto pair = models::makeAtmMatrix(m.asOf, m.curves, m.cube, {5}, {2, 10});
        const models::ExerciseRight& shortSwap = pair[0].right;
        const models::ExerciseRight& longSwap = pair[1].right;
        const double df = m.curves.discount().df(longSwap.time);
        const double curveSpread = longSwap.forwardRate() - shortSwap.forwardRate();
        Csv c(out, "cms_spread.csv",
              "model,rho,sigma_ratio,attainable,short_vol_error_bp,long_vol_error_bp,rate_correlation,spread_vol_bp,forward_spread_bp,price_bp,price_minus25_bp,price_plus25_bp");
        auto report2 = [&](const char* name, const G2pp& g) {
            auto err = [&](const models::CalibrationInstrument& ci) {
                return (models::impliedNormalVol(ci.right, models::g2ppSwaptionPrice(g, ci.right)) - ci.marketVol) * kBp;
            };
            const pricers::CmsSpreadResult atm = pricers::priceCmsSpread(g, longSwap, shortSwap, curveSpread, true, df);
            std::fprintf(c.f, "%s,%.2f,%.4f,1,%.1e,%.1e,%.4f,%.3f,%.3f,%.3f,%.3f,%.3f\n", name, g.rho(), g.sigma2() / g.sigma1(), err(pair[0]), err(pair[1]), atm.rateCorrelation,
                         atm.spreadVolBp, atm.forwardSpread * kBp, atm.price * kBp, pricers::priceCmsSpread(g, longSwap, shortSwap, curveSpread - 0.0025, true, df).price * kBp,
                         pricers::priceCmsSpread(g, longSwap, shortSwap, curveSpread + 0.0025, true, df).price * kBp);
        };
        report2("one factor", models::calibrateOneFactorToPair(pair[0], pair[1]));
        for (double rho : {-0.97, -0.95, -0.90, -0.85, -0.80, -0.75, -0.60, -0.30, 0.0}) {
            try {
                report2("two factors", models::calibrateG2ppToPair(fit2.model.kappa1(), fit2.model.kappa2(), rho, pair[0], pair[1], 0.5));
            } catch (const CalibrationError&) {
                std::fprintf(c.f, "two factors,%.2f,,0,,,,,,,,\n", rho);
            }
        }
        report2("fitted to the ATM matrix", fit2.model);
        Csv t(out, "cms_spread_info.csv", "expiry_years,df_expiry,curve_spread_bp,short_forward_pct,long_forward_pct,short_vol_bp,long_vol_bp,ms_per_price");
        std::fprintf(t.f, "%.4f,%.8f,%.3f,%.4f,%.4f,%.3f,%.3f,%.2f\n", longSwap.time, df, curveSpread * kBp, shortSwap.forwardRate() * 100.0, longSwap.forwardRate() * 100.0, pair[0].marketVol * kBp,
                     pair[1].marketVol * kBp, timeMs([&] { (void)pricers::priceCmsSpread(fit2.model, longSwap, shortSwap, curveSpread, true, df); }));
    }
}

// ---------------------------------------------------------------------------------------------------
// Exposure
// ---------------------------------------------------------------------------------------------------

void exposure(const std::filesystem::path& out, const models::UsdMarket& m) {
    Csv c(out, "exposure.csv", "time,ee_closed_form_bp,ee_monte_carlo_bp,mc_se_bp,pfe_975_bp,ene_bp");
    const instruments::Swap s = swap(m, 10, 0.035);
    const CalibratedBermudan b = pricers::calibrateBermudan(m, s, 1.0, 0, kKappa);
    for (const auto& p : pricers::swapExposureProfile(b.problem.model(), s, m.curves))
        std::fprintf(c.f, "%.3f,%.3f,%.3f,%.3f,%.3f,%.3f\n", p.time, p.expectedExposure * kBp, p.expectedExposureMc * kBp, p.mcStandardError * kBp,
                     p.potentialExposure * kBp, p.expectedNegativeExposure * kBp);
}

// ---------------------------------------------------------------------------------------------------
// Smile: which strike is the model calibrated at?
// ---------------------------------------------------------------------------------------------------

void smileImpact(const std::filesystem::path& out, const models::UsdMarket& m) {
    const auto withRate = [&](double k) { return swap(m, 10, k); };
    Csv c(out, "smile_bermudan.csv", "trade,calibration,bermudan_bp,best_european_bp,max_strip_error_bp");
    Csv d(out, "smile_strikes.csv", "trade,expiry,forward_pct,boundary_pct,vol_trade_strike_bp,vol_atm_bp,vol_boundary_bp");
    for (const auto& [name, strike, payer] : {std::tuple{"payer 3.50%", 0.035, true}, std::tuple{"payer 4.50% (OTM)", 0.045, true}, std::tuple{"payer 2.50% (ITM)", 0.025, true}}) {
        const instruments::Swap trade = swap(m, 10, strike, payer);
        const CalibratedBermudan base = pricers::calibrateBermudan(m, trade, 1.0, 0, kKappa);
        const std::size_t n = base.problem.size();
        auto price = [&](const CalibratedBermudan& b) { return pricers::priceBermudanPde(b.problem, pdeGrid(800)); };
        auto report = [&](const char* how, const CalibratedBermudan& b) {
            const std::vector<double> e = b.problem.europeanValues();
            std::fprintf(c.f, "%s,%s,%.3f,%.3f,%.1e\n", name, how, price(b).value * kBp, *std::max_element(e.begin(), e.end()) * kBp, b.calibration.maxAbsErrorBp);
        };
        report("trade strike", base);
        std::vector<double> atm;
        for (const auto& r : base.problem.rights()) atm.push_back(r.forwardRate());
        const CalibratedBermudan atTheMoney = pricers::calibrateBermudanAtStrikes(m, trade, 1.0, 0, kKappa, withRate, atm);
        report("at the money", atTheMoney);
        // The exercise boundary depends on the model, which depends on the calibration strikes: iterate.
        std::vector<double> boundary(n, strike);
        CalibratedBermudan atBoundary = base;
        for (int it = 0; it < 4; ++it) {
            const std::vector<double> rate = price(atBoundary).boundaryRate;
            for (std::size_t i = 0; i < n; ++i)
                if (std::isfinite(rate[i])) boundary[i] = rate[i];
            atBoundary = pricers::calibrateBermudanAtStrikes(m, trade, 1.0, 0, kKappa, withRate, boundary);
        }
        report("exercise boundary", atBoundary);
        for (std::size_t i = 0; i < n; ++i)
            std::fprintf(d.f, "%s,%.3f,%.4f,%.4f,%.3f,%.3f,%.3f\n", name, base.problem.rights()[i].time, atm[i] * 100.0, boundary[i] * 100.0, base.calibration.rows[i].marketVolBp,
                         atTheMoney.calibration.rows[i].marketVolBp, atBoundary.calibration.rows[i].marketVolBp);
    }
}

// ---------------------------------------------------------------------------------------------------
// Marking the mean reversion
// ---------------------------------------------------------------------------------------------------

void meanReversion(const std::filesystem::path& out, const models::UsdMarket& m) {
    const instruments::Swap trade = swap(m, 10, 0.035);
    std::vector<models::ExerciseRight> rights;
    for (core::Date d : models::coterminalExerciseDates(trade, 1.0, 0, "USD-SIFMA")) rights.push_back(models::makeExerciseRight(trade, m.curves, d));
    const auto strip = models::makeStrip(rights, m.cube);
    for (int tenor : {1, 2}) {
        const auto second = models::makeAtmMatrix(m.asOf, m.curves, m.cube, {1, 2, 3, 4, 5, 6, 7, 8, 9}, {tenor});
        const models::MeanReversionCalibration cal = models::calibrateLgmMeanReversion(strip, second);
        const std::string tag = std::to_string(tenor) + "y";
        {
            Csv c(out, ("kappa_scan_" + tag + ".csv").c_str(), "kappa,rms_bp");
            for (std::size_t i = 0; i < cal.scanKappa.size(); ++i) std::fprintf(c.f, "%.4f,%.4f\n", cal.scanKappa[i], cal.scanRmsBp[i]);
        }
        {
            Csv c(out, ("kappa_fit_" + tag + ".csv").c_str(), "label,market_vol_bp,error_bp");
            for (std::size_t i = 0; i < second.size(); ++i) std::fprintf(c.f, "%s,%.3f,%.3f\n", second[i].right.label.c_str(), second[i].marketVol * kBp, cal.secondaryErrorsBp[i]);
        }
        Csv c(out, ("kappa_result_" + tag + ".csv").c_str(), "kappa,rms_bp,at_bound,strip_error_bp,bermudan_bp,bermudan_at_reference_bp");
        const double marked = pricers::priceBermudanPde(BermudanProblem(cal.calibration.model, rights), pdeGrid(800)).value;
        const double reference = pricers::priceBermudanPde(pricers::calibrateBermudan(m, trade, 1.0, 0, kKappa).problem, pdeGrid(800)).value;
        std::fprintf(c.f, "%.5f,%.4f,%d,%.1e,%.3f,%.3f\n", cal.kappa, cal.secondaryRmsBp, cal.atBound ? 1 : 0, cal.calibration.maxAbsErrorBp, marked * kBp, reference * kBp);
    }
}

// ---------------------------------------------------------------------------------------------------
// XVA
// ---------------------------------------------------------------------------------------------------

void xvaSection(const std::filesystem::path& out, const models::UsdMarket& m) {
    using xva::CdsQuote;
    using xva::HazardCurve;
    const double recovery = 0.4, fundingSpread = 0.0060;
    const std::vector<CdsQuote> cptyQuotes{{1.0, 0.0080}, {3.0, 0.0105}, {5.0, 0.0125}, {7.0, 0.0140}, {10.0, 0.0150}};
    const std::vector<CdsQuote> ownQuotes{{1.0, 0.0040}, {3.0, 0.0055}, {5.0, 0.0065}, {7.0, 0.0072}, {10.0, 0.0078}};
    auto discount = [&](double t) { return m.curves.discount().df(t); };
    xva::XvaInputs in;
    in.counterparty = HazardCurve::bootstrap(cptyQuotes, recovery, discount);
    in.own = HazardCurve::bootstrap(ownQuotes, recovery, discount);
    in.counterpartyRecovery = in.ownRecovery = recovery;
    in.fundingSpread = fundingSpread;
    {
        Csv c(out, "xva_credit.csv", "name,maturity,spread_bp,hazard_pct,survival_pct,repriced_spread_bp");
        auto dump = [&](const char* name, const std::vector<CdsQuote>& q, const HazardCurve& h) {
            for (std::size_t i = 0; i < q.size(); ++i)
                std::fprintf(c.f, "%s,%.0f,%.1f,%.4f,%.3f,%.10f\n", name, q[i].maturityYears, q[i].spread * kBp, h.hazards()[i] * 100.0, h.survival(q[i].maturityYears) * 100.0,
                             xva::cdsParSpread(h, q[i].maturityYears, recovery, discount) * kBp);
        };
        dump("counterparty", cptyQuotes, in.counterparty);
        dump("own", ownQuotes, in.own);
    }

    const instruments::Swap ref = swap(m, 10, 0.035);
    const models::Lgm1F model = pricers::calibrateBermudan(m, ref, 1.0, 0, kKappa).problem.model();

    // --- One swap: closed form against the engine ----------------------------------------------------
    {
        const xva::SwapExposure exact = xva::swapExposureClosedForm(model, ref, m.curves);
        xva::ExposureOptions eo;
        eo.paths = 100000;
        const xva::ExposureProfile mc = xva::simulateExposure(model, {ref}, m.curves, exact.times, std::nullopt, eo);
        Csv c(out, "xva_swap.csv", "time,epe_closed_bp,epe_mc_bp,mc_se_bp,ene_closed_bp,ene_mc_bp");
        for (std::size_t i = 0; i < exact.times.size(); ++i)
            std::fprintf(c.f, "%.3f,%.3f,%.3f,%.3f,%.3f,%.3f\n", exact.times[i], exact.epe[i] * kBp, mc.epe[i] * kBp, mc.epeStandardError[i] * kBp, exact.ene[i] * kBp,
                         mc.ene[i] * kBp);
        const xva::XvaResult a = xva::computeXva(exact.times, exact.epe, exact.ene, in), b = xva::computeXva(mc.times, mc.epe, mc.ene, in);
        Csv s(out, "xva_swap_summary.csv", "method,cva_bp,dva_bp,fca_bp,fba_bp");
        std::fprintf(s.f, "closed form,%.4f,%.4f,%.4f,%.4f\n", a.cva * kBp, a.dva * kBp, a.fca * kBp, a.fba * kBp);
        std::fprintf(s.f, "Monte Carlo,%.4f,%.4f,%.4f,%.4f\n", b.cva * kBp, b.dva * kBp, b.fca * kBp, b.fba * kBp);
    }

    // --- A netting set, with and without collateral ---------------------------------------------------
    const struct {
        int years;
        double rate;
        bool payer;
        double notional;
    } book[] = {{10, 0.0350, true, 1.0}, {7, 0.0340, false, 0.8}, {5, 0.0360, true, 1.2}, {3, 0.0330, false, 1.5},
                {2, 0.0350, true, 1.0},  {10, 0.0370, false, 0.6}, {7, 0.0320, true, 0.5}, {5, 0.0355, false, 0.9}};
    std::vector<instruments::Swap> trades;
    {
        Csv c(out, "xva_book.csv", "trade,tenor_years,fixed_rate_pct,side,notional,value_bp");
        int k = 0;
        for (const auto& b : book) {
            trades.push_back(swap(m, b.years, b.rate, b.payer, b.notional));
            const xva::SwapExposure e = xva::swapExposureClosedForm(model, trades.back(), m.curves);
            std::fprintf(c.f, "%d,%d,%.2f,%s,%.1f,%.2f\n", ++k, b.years, b.rate * 100.0, b.payer ? "pay fixed" : "receive fixed", b.notional, (e.epe[0] + e.ene[0]) * kBp);
        }
    }
    const std::vector<double> times = xva::exposureTimes(trades, m.curves, 0.25);
    xva::ExposureOptions eo;
    eo.paths = 50000;
    const double tenDays = 10.0 / 252.0;
    const struct {
        const char* name;
        xva::Csa csa;
    } agreements[] = {{"zero threshold; 10-day margin period", {0.0, 0.0, tenDays}},
                      {"threshold 50 bp of notional", {0.0050, 0.0, tenDays}},
                      {"minimum transfer 10 bp", {0.0, 0.0010, tenDays}},
                      {"zero threshold; 20-day margin period", {0.0, 0.0, 2.0 * tenDays}}};
    std::vector<xva::ExposureProfile> runs;
    double engineMs = 0.0;
    for (const auto& a : agreements) {
        const auto t0 = std::chrono::steady_clock::now();
        runs.push_back(xva::simulateExposure(model, trades, m.curves, times, a.csa, eo));
        engineMs = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
    }
    const xva::ExposureProfile& base = runs[0];
    {
        Csv c(out, "xva_portfolio.csv", "time,epe_bp,ene_bp,pfe_bp,gross_epe_bp,csa_zero_bp,csa_threshold_bp,csa_mta_bp,csa_20d_bp,coupon_date");
        // A coupon date: some trade's fixed period starts here (its previous coupon has just rolled off).
        std::vector<double> starts;
        for (const auto& t : trades)
            for (double s : xva::periodStartTimes(t, m.curves)) starts.push_back(s);
        for (std::size_t i = 0; i < times.size(); ++i) {
            bool coupon = false;
            for (double s : starts)
                if (s > 0.1 && std::fabs(s - times[i]) < 1e-6) coupon = true;
            std::fprintf(c.f, "%.4f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%d\n", times[i], base.epe[i] * kBp, base.ene[i] * kBp, base.pfe[i] * kBp, base.grossEpe[i] * kBp,
                         runs[0].collateralisedEpe[i] * kBp, runs[1].collateralisedEpe[i] * kBp, runs[2].collateralisedEpe[i] * kBp, runs[3].collateralisedEpe[i] * kBp,
                         coupon ? 1 : 0);
        }
    }
    {
        Csv c(out, "xva_table.csv", "case,cva_bp,dva_bp,fca_bp,fba_bp");
        auto row = [&](const char* name, const std::vector<double>& epe, const std::vector<double>& ene) {
            const xva::XvaResult r = xva::computeXva(times, epe, ene, in);
            std::fprintf(c.f, "%s,%.4f,%.4f,%.4f,%.4f\n", name, r.cva * kBp, r.dva * kBp, r.fca * kBp, r.fba * kBp);
        };
        // With no netting agreement each trade is its own claim: CVA on the sum of the trades' exposures.
        std::vector<double> none(times.size(), 0.0);
        const xva::XvaResult gross = xva::computeXva(times, base.grossEpe, none, in);
        std::fprintf(c.f, "no netting agreement,%.4f,,%.4f,\n", gross.cva * kBp, gross.fca * kBp);
        row("netted; no collateral", base.epe, base.ene);
        for (std::size_t k = 0; k < runs.size(); ++k) row(agreements[k].name, runs[k].collateralisedEpe, runs[k].collateralisedEne);
        Csv t(out, "xva_engine.csv", "trades,exposure_dates,paths,ms");
        std::fprintf(t.f, "%zu,%zu,%zu,%.0f\n", trades.size(), times.size(), eo.paths, engineMs);
    }

    // --- The Bermudan: an option that turns into a swap ------------------------------------------------
    {
        const CalibratedBermudan b = pricers::calibrateBermudan(m, ref, 1.0, 0, kKappa);
        xva::BermudanExposureOptions bo;
        bo.paths = 200000;
        const xva::BermudanExposure e = xva::bermudanExposure(b.problem, bo);
        const xva::SwapExposure s = xva::swapExposureClosedForm(model, ref, m.curves);
        const std::vector<double> europeans = b.problem.europeanValues();
        Csv c(out, "xva_bermudan.csv", "time,epe_bp,se_bp,option_epe_bp,swap_epe_bp,ene_bp,exercised_pct,underlying_swap_epe_bp,european_bp");
        for (std::size_t i = 0; i < e.times.size(); ++i)
            std::fprintf(c.f, "%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.2f,%.3f,%.3f\n", e.times[i], e.epe[i] * kBp, e.epeStandardError[i] * kBp, e.optionEpe[i] * kBp,
                         (e.epe[i] - e.optionEpe[i]) * kBp, e.ene[i] * kBp, e.exercisedBy[i] * 100.0, (i == 0 ? s.epe[0] : s.epe[i + 1]) * kBp,
                         i == 0 ? 0.0 : europeans[i - 1] * kBp);
        // The underlying swap's own profile on the same dates (dropping its spot date) for a like-for-like CVA.
        std::vector<double> swapEpe{s.epe[0]}, swapEne{s.ene[0]};
        for (std::size_t i = 2; i < s.times.size(); ++i) {
            swapEpe.push_back(s.epe[i]);
            swapEne.push_back(s.ene[i]);
        }
        const xva::XvaResult rb = xva::computeXva(e.times, e.epe, e.ene, in), rs = xva::computeXva(e.times, swapEpe, swapEne, in);
        Csv t(out, "xva_bermudan_summary.csv", "grid_price_bp,simulated_price_bp,simulated_se_bp,pde_price_bp,cva_bp,dva_bp,fca_bp,swap_cva_bp,swap_dva_bp");
        std::fprintf(t.f, "%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n", e.price * kBp, e.simulatedPrice * kBp, e.simulatedPriceStandardError * kBp,
                     pricers::priceBermudanPde(b.problem, pdeGrid(800)).value * kBp, rb.cva * kBp, rb.dva * kBp, rb.fca * kBp, rs.cva * kBp, rs.dva * kBp);
    }

    // --- Sensitivities ---------------------------------------------------------------------------------
    {
        const xva::SwapExposure e = xva::swapExposureClosedForm(model, ref, m.curves);
        auto swapCva = [&](const HazardCurve& h) { return xva::creditAdjustment(e.times, e.epe, h, recovery); };
        auto bookCva = [&](const HazardCurve& h) { return xva::creditAdjustment(times, base.epe, h, recovery); };
        const auto a = xva::cs01Ladder(cptyQuotes, recovery, discount, swapCva), b = xva::cs01Ladder(cptyQuotes, recovery, discount, bookCva);
        Csv c(out, "xva_cs01.csv", "quote,swap_cva_cs01_bp,book_cva_cs01_bp");
        for (std::size_t i = 0; i < a.size(); ++i) std::fprintf(c.f, "%s,%.5f,%.5f\n", a[i].id.c_str(), a[i].value * kBp, b[i].value * kBp);
    }
    {
        const HazardCurve& hazard = in.counterparty;
        auto frozen = [&](const models::UsdMarket& mk) { return xva::swapCvaClosedForm(model, ref, mk.curves, hazard, recovery); };
        auto recalibrated = [&](const models::UsdMarket& mk) {
            return xva::swapCvaClosedForm(pricers::calibrateBermudan(mk, ref, 1.0, 0, kKappa).problem.model(), ref, mk.curves, hazard, recovery);
        };
        const auto adj = xva::swapCvaDeltaAdjoint(m, model, ref, hazard, recovery);
        const auto bumpFrozen = pricers::parDeltaBump(m, frozen), bumpRecal = pricers::parDeltaBump(m, recalibrated);
        Csv c(out, "xva_delta.csv", "quote,adjoint,bump_model_fixed,bump_recalibrated");
        for (std::size_t i = 0; i < adj.size(); ++i)
            std::fprintf(c.f, "%s,%.6f,%.6f,%.6f\n", adj[i].id.c_str(), adj[i].value * kBp, bumpFrozen[i].value * kBp, bumpRecal[i].value * kBp);
        Csv t(out, "xva_delta_timing.csv", "method,ms");
        std::fprintf(t.f, "one CVA,%.4f\n", timeMs([&] { (void)frozen(m); }));
        std::fprintf(t.f, "adjoint,%.4f\n", timeMs([&] { (void)xva::swapCvaDeltaAdjoint(m, model, ref, hazard, recovery); }));
        std::fprintf(t.f, "bump (model fixed),%.4f\n", timeMs([&] { (void)pricers::parDeltaBump(m, frozen); }));
        std::fprintf(t.f, "bump (recalibrated),%.4f\n", timeMs([&] { (void)pricers::parDeltaBump(m, recalibrated); }));
    }

    // --- Wrong-way risk: default that depends on rates ---------------------------------------------------
    {
        xva::ExposureOptions wo;
        wo.paths = 50000;
        wo.keepPaths = true;
        const instruments::Swap payer = swap(m, 10, 0.035, true), receiver = swap(m, 10, 0.035, false);
        const std::vector<double> swapTimes = xva::exposureTimes({payer}, m.curves, 0.25);
        const xva::ExposureProfile pp = xva::simulateExposure(model, {payer}, m.curves, swapTimes, std::nullopt, wo);
        const xva::ExposureProfile pr = xva::simulateExposure(model, {receiver}, m.curves, swapTimes, std::nullopt, wo);
        const xva::ExposureProfile pb = xva::simulateExposure(model, trades, m.curves, times, std::nullopt, wo);
        Csv c(out, "xva_wwr.csv", "beta,hazard_multiple_per_sd,payer_cva_bp,payer_se_bp,receiver_cva_bp,receiver_se_bp,book_cva_bp,book_se_bp,worst_survival_error");
        for (double beta : {-1.0, -0.5, -0.25, 0.0, 0.25, 0.5, 1.0}) {
            const xva::WrongWayResult a = xva::wrongWayCva(pp, model, in.counterparty, recovery, beta), b = xva::wrongWayCva(pr, model, in.counterparty, recovery, beta),
                                      d = xva::wrongWayCva(pb, model, in.counterparty, recovery, beta);
            double worst = 0.0;
            for (const auto* r : {&a, &b, &d})
                for (double e : r->survivalError) worst = std::max(worst, std::fabs(e));
            std::fprintf(c.f, "%.2f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.1e\n", beta, std::exp(beta), a.cva * kBp, a.standardError * kBp, b.cva * kBp, b.standardError * kBp,
                         d.cva * kBp, d.standardError * kBp, worst);
        }
        const xva::WrongWayResult up = xva::wrongWayCva(pp, model, in.counterparty, recovery, 1.0), down = xva::wrongWayCva(pp, model, in.counterparty, recovery, -1.0);
        Csv p(out, "xva_wwr_profile.csv", "time,epe_bp,epe_given_default_wrong_way_bp,epe_given_default_right_way_bp");
        for (std::size_t i = 1; i < up.times.size(); ++i)
            std::fprintf(p.f, "%.4f,%.3f,%.3f,%.3f\n", up.times[i], up.epe[i] * kBp, up.epeGivenDefault[i] * kBp, down.epeGivenDefault[i] * kBp);
    }

    // --- Two factors: the same book and the same Bermudan when rates can decorrelate --------------------
    {
        using models::G2pp;
        // Shape from the ATM matrix (section 9.2), level from the 10-year co-terminal strip, as for the
        // one-factor model. Every model here reprices that strip, so the 10-year swap's own exposure on
        // its coupon dates is the same in all of them.
        models::G2ppFitOptions fo;
        fo.kappa1 = 0.6;
        fo.rho = -0.9;
        const G2pp shape = models::fitG2pp(models::makeAtmMatrix(m.asOf, m.curves, m.cube, {1, 2, 3, 5, 7, 10}, {1, 2, 5, 10, 20}), fo).model;
        const std::vector<double> couponTimes = xva::periodStartTimes(ref, m.curves);
        Csv c(out, "xva_twofactor.csv", "model,rho,swap_cva_bp,book_cva_bp,book_gross_cva_bp,book_dva_bp");
        xva::ExposureOptions e2;
        e2.paths = 50000;
        auto cvaOf = [&](const std::vector<double>& t, const std::vector<double>& epe) { return xva::creditAdjustment(t, epe, in.counterparty, recovery) * kBp; };
        {
            const xva::ExposureProfile a = xva::simulateExposure(model, {ref}, m.curves, couponTimes, std::nullopt, e2);
            std::fprintf(c.f, "one factor,,%.4f,%.4f,%.4f,%.4f\n", cvaOf(a.times, a.epe), cvaOf(times, base.epe), cvaOf(times, base.grossEpe),
                         xva::computeXva(times, base.epe, base.ene, in).dva * kBp);
        }
        auto row = [&](const char* name, double rho) {
            const pricers::CalibratedBermudanG2pp cal =
                pricers::calibrateBermudanG2pp(m, ref, 1.0, 0, G2pp(shape.kappa1(), shape.kappa2(), shape.sigma1(), shape.sigma2(), rho));
            const G2pp& g = cal.calibration.model;
            const xva::ExposureProfile a = xva::simulateExposure(g, {ref}, m.curves, couponTimes, std::nullopt, e2);
            const xva::ExposureProfile b = xva::simulateExposure(g, trades, m.curves, times, std::nullopt, e2);
            std::fprintf(c.f, "%s,%.2f,%.4f,%.4f,%.4f,%.4f\n", name, rho, cvaOf(a.times, a.epe), cvaOf(b.times, b.epe), cvaOf(b.times, b.grossEpe),
                         xva::computeXva(b.times, b.epe, b.ene, in).dva * kBp);
            return cal;
        };
        for (double rho : {0.5, 0.0, -0.5, -0.9}) row("two factors", rho);
        const pricers::CalibratedBermudanG2pp fitted = row("two factors; fitted shape", shape.rho());

        xva::BermudanExposureG2ppOptions bo;
        bo.trainingPaths = 100000;
        bo.paths = 400000;
        const xva::BermudanExposure e = xva::bermudanExposureG2pp(fitted.calibration.model, fitted.rights, bo);
        Csv b(out, "xva_bermudan_2f.csv", "time,epe_bp,se_bp,option_epe_bp,swap_epe_bp,ene_bp");
        for (std::size_t i = 0; i < e.times.size(); ++i)
            std::fprintf(b.f, "%.3f,%.3f,%.3f,%.3f,%.3f,%.3f\n", e.times[i], e.epe[i] * kBp, e.epeStandardError[i] * kBp, e.optionEpe[i] * kBp, (e.epe[i] - e.optionEpe[i]) * kBp,
                         e.ene[i] * kBp);
        const xva::XvaResult r = xva::computeXva(e.times, e.epe, e.ene, in);
        Csv t(out, "xva_bermudan_2f_summary.csv", "price_bp,price_se_bp,cva_bp,dva_bp,rho,kappa1,kappa2");
        std::fprintf(t.f, "%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n", e.price * kBp, e.simulatedPriceStandardError * kBp, r.cva * kBp, r.dva * kBp, shape.rho(), shape.kappa1(),
                     shape.kappa2());
    }
}

} // namespace

int main(int argc, char** argv) {
    if (argc < 2) {
        std::fprintf(stderr, "usage: rl_analysis <output-dir> [section...]\n  sections: curves swaprisk vol lgm bermudan exposure twofactor smile kappa xva (default: all)\n");
        return 2;
    }
    try {
        const std::filesystem::path out = argv[1];
        std::filesystem::create_directories(out);
        auto wanted = [&](const char* name) {
            if (argc == 2) return true;
            for (int i = 2; i < argc; ++i)
                if (std::string(argv[i]) == name) return true;
            return false;
        };
        const models::UsdMarket m = models::loadUsdMarket();
        if (wanted("curves")) curves(out);
        if (wanted("swaprisk")) swapRisk(out, m);
        if (wanted("vol")) volCube(out, m);
        if (wanted("lgm")) lgm(out, m);
        if (wanted("bermudan")) bermudan(out, m);
        if (wanted("exposure")) exposure(out, m);
        if (wanted("twofactor")) twoFactor(out, m);
        if (wanted("smile")) smileImpact(out, m);
        if (wanted("kappa")) meanReversion(out, m);
        if (wanted("xva")) xvaSection(out, m);
        return 0;
    } catch (const std::exception& e) {
        std::fprintf(stderr, "rl_analysis: %s\n", e.what());
        return 1;
    }
}
