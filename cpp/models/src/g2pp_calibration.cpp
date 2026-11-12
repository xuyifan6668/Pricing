#include "rl/models/g2pp_calibration.hpp"

#include "rl/core/error.hpp"
#include "rl/core/tenor.hpp"
#include "rl/instruments/builder.hpp"
#include "rl/math/rootfind.hpp"
#include "rl/models/g2pp_swaption.hpp"
#include "rl/models/lgm_swaption.hpp"
#include "rl/vol/cube.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <variant>

namespace rl::models {

std::vector<CalibrationInstrument> makeAtmMatrix(core::Date asOf, const market::CurveSet<double>& curves, const market::SwaptionCube<double>& cube,
                                                 const std::vector<int>& expiryYears, const std::vector<int>& tenorYears) {
    std::vector<CalibrationInstrument> out;
    for (int e : expiryYears) {
        for (int n : tenorYears) {
            RL_REQUIRE(e > 0 && n > 0, "makeAtmMatrix: expiries and tenors must be positive whole years");
            // The swap [e, e + n] is the tail of a spot-starting (e + n)-year swap: exercise at the start
            // of its fixed period number e. Built twice, the second time at its own forward rate.
            auto right = [&](double rate) {
                const instruments::Swap swap =
                    instruments::OisSwapBuilder::usdSofr(asOf).id("atm").tenor(core::Tenor::years(e + n)).fixedRate(rate).notional(1.0).payFixed(true).build();
                const auto& cfs = swap.fixedLeg().cashflows();
                RL_REQUIRE(static_cast<std::size_t>(e) < cfs.size(), "makeAtmMatrix: the fixed leg has fewer periods than the expiry");
                return makeExerciseRight(swap, curves, std::get<instruments::FixedCoupon>(cfs[static_cast<std::size_t>(e)]).accrualStart);
            };
            CalibrationInstrument ci;
            ci.right = right(right(0.03).forwardRate());
            ci.right.label = std::to_string(e) + "Yx" + std::to_string(n) + "Y";
            ci.marketVol = vol::cubeAtmNormalVol<double>(cube, ci.right.time, ci.right.tenorYears);
            ci.marketPrice = marketSwaptionPrice(ci.right, ci.marketVol);
            out.push_back(std::move(ci));
        }
    }
    return out;
}

namespace {

const double kTwoPi = 6.283185307179586;

/// Model normal vol in bp. At the money the Bachelier price is A σ √(t/2π), so the vol is read off
/// directly; the same expression is the fallback where the inverter has no bracket (the optimiser
/// can visit parameters where a swaption is worth its intrinsic value, or absurdly much).
double modelVolBp(const G2pp& model, const CalibrationInstrument& ci) {
    const ExerciseRight& r = ci.right;
    const double price = g2ppSwaptionPrice(model, r);
    const double timeValue = std::max(price - std::max(r.forwardValue(), 0.0), 0.0);
    const double atm = timeValue / r.annuityValue() * std::sqrt(kTwoPi / r.time) * 1e4;
    if (std::fabs(r.forwardRate() - r.strike) < 1e-10) return atm;
    try {
        return impliedNormalVol(r, price) * 1e4;
    } catch (const Error&) {
        return atm;
    }
}

/// Parameter layout: [κ_1, (κ_2, log(σ_2/σ_1), ρ), log s_1 … log s_m]. σ_1 is fixed at 1% and the
/// level is carried by the scale; with no scale nodes the single "scale" is the constant level.
struct Layout {
    int factors;
    std::vector<double> scaleTimes; ///< empty: constant parameters
    int shape() const { return factors == 1 ? 1 : 4; }
    int levels() const { return scaleTimes.empty() ? 1 : static_cast<int>(scaleTimes.size()); }
    int size() const { return shape() + levels(); }

    G2pp model(const math::Vector& x) const {
        const double sigma1 = 0.01;
        const double kappa2 = factors == 1 ? 1.0 : x[1];
        const double sigma2 = factors == 1 ? 0.0 : sigma1 * std::exp(x[2]);
        const double rho = factors == 1 ? 0.0 : x[3];
        if (scaleTimes.empty()) {
            const double level = std::exp(x[shape()]);
            return G2pp(x[0], kappa2, sigma1 * level, sigma2 * level, rho);
        }
        std::vector<double> scales(scaleTimes.size());
        for (std::size_t j = 0; j < scales.size(); ++j) scales[j] = std::exp(x[shape() + static_cast<int>(j)]);
        return G2pp(x[0], kappa2, sigma1, sigma2, rho, scaleTimes, scales);
    }
};

} // namespace

G2ppFit fitG2pp(const std::vector<CalibrationInstrument>& instruments, const G2ppFitOptions& opt) {
    RL_REQUIRE(opt.factors == 1 || opt.factors == 2, "fitG2pp: factors must be 1 or 2");
    Layout layout{opt.factors, {}};
    if (opt.timeDependentScale) {
        for (const CalibrationInstrument& ci : instruments) layout.scaleTimes.push_back(ci.right.time);
        std::sort(layout.scaleTimes.begin(), layout.scaleTimes.end());
        // Expiries of the same grid row differ by a day or two of date rolling: one node per row.
        std::vector<double> distinct;
        for (double t : layout.scaleTimes)
            if (distinct.empty() || t > distinct.back() + 0.05) distinct.push_back(t);
            else distinct.back() = t;
        layout.scaleTimes = distinct;
    }
    const int n = layout.size(), ns = layout.shape();
    RL_REQUIRE(instruments.size() > static_cast<std::size_t>(n), "fitG2pp: need more instruments than parameters");

    math::Vector x0(n), lo(n), hi(n);
    x0[0] = opt.kappa1;
    lo[0] = opt.kappaMin;
    hi[0] = opt.kappaMax;
    if (opt.factors == 2) {
        x0[1] = opt.kappa2;
        lo[1] = opt.kappaMin;
        hi[1] = opt.kappaMax;
        x0[2] = std::log(opt.sigma2 / opt.sigma1);
        lo[2] = std::log(1e-3);
        hi[2] = std::log(1e3);
        x0[3] = opt.rho;
        lo[3] = -opt.rhoBound;
        hi[3] = opt.rhoBound;
    }
    for (int j = ns; j < n; ++j) {
        x0[j] = std::log(opt.sigma1 / 0.01);
        lo[j] = std::log(1e-3);
        hi[j] = std::log(1e2);
    }
    auto residuals = [&](const math::Vector& x) {
        const G2pp model = layout.model(x);
        math::Vector r(static_cast<Eigen::Index>(instruments.size()));
        for (std::size_t i = 0; i < instruments.size(); ++i) r[static_cast<Eigen::Index>(i)] = modelVolBp(model, instruments[i]) - instruments[i].marketVol * 1e4;
        return r;
    };
    math::LmOptions lm;
    lm.lowerBounds = lo;
    lm.upperBounds = hi;
    lm.maxIter = opt.maxIterations;
    const math::LmResult res = math::levenbergMarquardt(residuals, {}, x0, lm);

    G2ppFit fit{layout.model(res.x), {}, 0.0, 0.0, res.iterations, res.status, {}};
    double sum = 0.0;
    for (std::size_t i = 0; i < instruments.size(); ++i) {
        const ExerciseRight& r = instruments[i].right;
        G2ppFitRow row{r.label, r.time, r.tenorYears, instruments[i].marketVol * 1e4, 0.0, 0.0};
        row.errorBp = res.residuals[static_cast<Eigen::Index>(i)];
        row.modelVolBp = row.marketVolBp + row.errorBp;
        sum += row.errorBp * row.errorBp;
        fit.maxAbsErrorBp = std::max(fit.maxAbsErrorBp, std::fabs(row.errorBp));
        fit.rows.push_back(std::move(row));
    }
    fit.rmsErrorBp = std::sqrt(sum / static_cast<double>(instruments.size()));
    const char* names[] = {"kappa1", "kappa2", "sigma2/sigma1", "rho"};
    for (int j = 0; j < n; ++j) {
        if (res.x[j] > lo[j] + 1e-9 * (1.0 + std::fabs(lo[j])) && res.x[j] < hi[j] - 1e-9 * (1.0 + std::fabs(hi[j]))) continue;
        fit.activeBounds.emplace_back(j < ns ? names[j] : "scale " + std::to_string(j - ns + 1));
    }
    return fit;
}

G2ppStripCalibration calibrateG2ppToStrip(const G2pp& shape, const std::vector<CalibrationInstrument>& strip) {
    RL_REQUIRE(!strip.empty(), "calibrateG2ppToStrip: empty strip");
    for (std::size_t i = 1; i < strip.size(); ++i)
        RL_REQUIRE(strip[i].right.time > strip[i - 1].right.time, "calibrateG2ppToStrip: the strip must be in increasing expiry order");

    std::vector<double> times, scales;
    std::vector<bool> clamped;
    for (const CalibrationInstrument& ci : strip) {
        times.push_back(ci.right.time);
        scales.push_back(1.0);
        // The price depends on the scales up to this expiry only, and rises with the last one.
        auto error = [&](double logScale) {
            scales.back() = std::exp(logScale);
            return g2ppSwaptionPrice(shape.withScale(times, scales), ci.right) - ci.marketPrice;
        };
        const double lo = std::log(1e-4);
        bool clamp = false;
        if (error(lo) >= 0.0) {
            clamp = scales.size() > 1; // the earlier intervals already give too much variance
            scales.back() = std::exp(lo);
        } else {
            // Widen the bracket from a scale of 2: the shape's own volatilities are of the right order,
            // and at absurd scales the price saturates at the value of the receive side.
            double hi = std::log(2.0);
            while (error(hi) <= 0.0) {
                hi += std::log(2.0);
                RL_REQUIRE(hi < std::log(200.0), "calibrateG2ppToStrip: market price of " + ci.right.label + " is not attainable");
            }
            math::RootOptions ro;
            ro.xTol = 1e-13;
            ro.maxIter = 200;
            scales.back() = std::exp(math::brent(error, lo, hi, ro).x);
        }
        clamped.push_back(clamp);
    }

    G2ppStripCalibration out{shape.withScale(times, scales), {}, 0.0, true};
    for (std::size_t i = 0; i < strip.size(); ++i) {
        const ExerciseRight& r = strip[i].right;
        G2ppStripRow row{r.label, r.time, r.tenorYears, strip[i].marketVol * 1e4, modelVolBp(out.model, strip[i]), 0.0, scales[i], clamped[i]};
        row.errorBp = row.modelVolBp - row.marketVolBp;
        out.maxAbsErrorBp = std::max(out.maxAbsErrorBp, std::fabs(row.errorBp));
        if (row.clamped) out.ok = false;
        out.rows.push_back(std::move(row));
    }
    if (out.maxAbsErrorBp > 0.01) out.ok = false;
    return out;
}

namespace {

/// The volatility level that reprices `ci` for a model shape given at unit level; the price rises with it.
double levelFor(const std::function<G2pp(double)>& atLevel, const CalibrationInstrument& ci) {
    auto error = [&](double logLevel) { return g2ppSwaptionPrice(atLevel(std::exp(logLevel)), ci.right) - ci.marketPrice; };
    // Walk the bracket up from a negligible level: at absurd volatilities the price saturates or overflows.
    const double lo = std::log(1e-4);
    RL_REQUIRE(error(lo) < 0.0, "G2++ pair calibration: market price of " + ci.right.label + " is at or below intrinsic");
    double hi = lo + std::log(2.0);
    while (!(error(hi) > 0.0)) {
        hi += std::log(2.0);
        RL_REQUIRE(hi < std::log(1e3), "G2++ pair calibration: market price of " + ci.right.label + " is not attainable");
    }
    math::RootOptions ro;
    ro.xTol = 1e-13;
    ro.maxIter = 200;
    return std::exp(math::brent(error, lo, hi, ro).x);
}

/// Solve a one-parameter family of shapes for the member that, once its level reprices `first`,
/// also reprices `second`. The second vol need not be monotone in the shape parameter, so there can
/// be more than one such member: the one closest to `preferred` is returned.
G2pp solvePair(const std::function<G2pp(double, double)>& build, double lo, double hi, double preferred, bool logScale, const CalibrationInstrument& first,
               const CalibrationInstrument& second, const char* who) {
    auto at = [&](double p) {
        const double shape = logScale ? std::exp(p) : p;
        return build(shape, levelFor([&](double level) { return build(shape, level); }, first));
    };
    auto error = [&](double p) { return modelVolBp(at(p), second) - second.marketVol * 1e4; };
    const double a = logScale ? std::log(lo) : lo, b = logScale ? std::log(hi) : hi, target = logScale ? std::log(preferred) : preferred;
    const int n = 60;
    double prev = a, ePrev = error(a), nearest = std::fabs(ePrev), root = 0.0;
    bool found = false;
    for (int i = 1; i <= n; ++i) {
        const double p = a + (b - a) * i / n, e = error(p);
        nearest = std::min(nearest, std::fabs(e));
        if ((ePrev < 0.0) != (e < 0.0)) {
            math::RootOptions ro;
            ro.xTol = 1e-12;
            ro.maxIter = 200;
            const double r = math::brent(error, prev, p, ro).x;
            if (!found || std::fabs(r - target) < std::fabs(root - target)) root = r;
            found = true;
        }
        prev = p;
        ePrev = e;
    }
    if (!found) throw CalibrationError(std::string(who) + ": no model in the range reprices both " + first.right.label + " and " + second.right.label, {0.0, nearest});
    return at(root);
}

} // namespace

G2pp calibrateG2ppToPair(double kappa1, double kappa2, double rho, const CalibrationInstrument& first, const CalibrationInstrument& second,
                         double preferredRatio, double ratioMin, double ratioMax) {
    RL_REQUIRE(ratioMin > 0.0 && ratioMax > ratioMin && preferredRatio > 0.0, "calibrateG2ppToPair: need 0 < ratioMin < ratioMax and a positive preferred ratio");
    return solvePair([&](double ratio, double level) { return G2pp(kappa1, kappa2, 0.01 * level, 0.01 * level * ratio, rho); }, ratioMin, ratioMax, preferredRatio,
                     true, first, second, "calibrateG2ppToPair");
}

G2pp calibrateOneFactorToPair(const CalibrationInstrument& first, const CalibrationInstrument& second, double kappaMin, double kappaMax) {
    RL_REQUIRE(kappaMax > kappaMin, "calibrateOneFactorToPair: empty mean-reversion range");
    return solvePair([&](double kappa, double level) { return G2pp(kappa, 1.0, 0.01 * level, 0.0, 0.0); }, kappaMin, kappaMax, 0.0, false, first, second,
                     "calibrateOneFactorToPair");
}

} // namespace rl::models
