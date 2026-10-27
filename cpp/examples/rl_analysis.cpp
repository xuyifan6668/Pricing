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
#include "rl/models/lgm_swaption.hpp"
#include "rl/math/rng.hpp"
#include "rl/models/usd_market.hpp"
#include "rl/pricers/bermudan_mc.hpp"
#include "rl/pricers/bermudan_pde.hpp"
#include "rl/pricers/exposure.hpp"
#include "rl/pricers/risk.hpp"
#include "rl/vol/arbitrage.hpp"
#include "rl/vol/sabr.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <functional>
#include <iostream>
#include <string>
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

} // namespace

int main(int argc, char** argv) {
    if (argc < 2) {
        std::fprintf(stderr, "usage: rl_analysis <output-dir> [section...]\n  sections: curves swaprisk vol lgm bermudan exposure (default: all)\n");
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
        return 0;
    } catch (const std::exception& e) {
        std::fprintf(stderr, "rl_analysis: %s\n", e.what());
        return 1;
    }
}
