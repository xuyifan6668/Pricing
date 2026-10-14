// implied normal and shifted-lognormal vols (double only). Sources and the bracketing argument are
// in bachelier.hpp and black.hpp; the root finders are math::newtonSafeguarded (Numerical Recipes 3e §9.4) and
// math::brent (§9.3).
#include "rl/vol/bachelier.hpp"
#include "rl/vol/black.hpp"

#include "rl/core/error.hpp"
#include "rl/math/normal.hpp"
#include "rl/math/rootfind.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>

namespace rl::vol {

namespace {
constexpr double kSqrt2Pi = 2.50662827463100050242;
constexpr double kTiny = 1e-300;

/// Time value = price − intrinsic, with a relative tolerance for a price that rounds a hair below intrinsic.
double timeValue(const char* who, double price, double intrinsic) {
    const double tv = price - intrinsic;
    const double tol = 1e-13 * std::max({std::fabs(price), std::fabs(intrinsic), 1e-12});
    RL_REQUIRE(tv >= -tol, std::string(who) + ": price " + std::to_string(price) + " is below intrinsic value " + std::to_string(intrinsic));
    return tv;
}
} // namespace

double bachelierImpliedVol(double price, double forward, double strike, double expiry, OptionType type, const ImpliedVolOptions& options) {
    RL_REQUIRE(std::isfinite(price) && std::isfinite(forward) && std::isfinite(strike), "bachelierImpliedVol: inputs must be finite");
    RL_REQUIRE(std::isfinite(expiry) && expiry > 0.0, "bachelierImpliedVol: expiry must be positive");
    const double m = forward - strike;
    const double intrinsic = std::max(optionSign(type) * m, 0.0);
    const double tv = timeValue("bachelierImpliedVol", price, intrinsic);
    if (tv <= 0.0) return 0.0;
    const double sqrtT = std::sqrt(expiry);
    const double x = std::fabs(m);
    if (x == 0.0) return tv * kSqrt2Pi / sqrtT; // ATM: C = s/√(2π) exactly

    // h(s) = s φ(x/s) − x Φ(−x/s): the OTM price as a function of s = σ√T; increasing, h'(s) = φ(x/s).
    // Bracket from s/√(2π) − x/2 ≤ h(s) ≤ s/√(2π).
    const double lo = tv * kSqrt2Pi;
    const double hi = (tv + 0.5 * x) * kSqrt2Pi;
    auto fdf = [&](double s) -> std::pair<double, double> {
        const double d = x / s;
        const double h = s * math::normalPdf(d) - x * math::normalCdf(-d);
        if (h <= kTiny) return {-800.0 - std::log(tv), x * x / (s * s * s)}; // underflow: log h ≈ −x²/(2s²)
        return {std::log(h / tv), math::normalPdf(d) / h};
    };
    math::RootOptions ro;
    ro.xTol = options.relTol * hi;
    ro.maxIter = options.maxIter;
    const double s0 = std::clamp(std::sqrt(lo * hi), lo, hi);
    if (lo >= hi) return lo / sqrtT;
    const math::RootResult r = math::newtonSafeguarded(fdf, s0, lo, hi, ro);
    return r.x / sqrtT;
}

double blackImpliedVol(double price, double forward, double strike, double expiry, double shift, OptionType type, const ImpliedVolOptions& options) {
    RL_REQUIRE(std::isfinite(price) && std::isfinite(forward) && std::isfinite(strike) && std::isfinite(shift), "blackImpliedVol: inputs must be finite");
    RL_REQUIRE(std::isfinite(expiry) && expiry > 0.0, "blackImpliedVol: expiry must be positive");
    const double f = forward + shift;
    const double k = strike + shift;
    RL_REQUIRE(f > 0.0, "blackImpliedVol: forward + shift must be positive");
    RL_REQUIRE(k > 0.0, "blackImpliedVol: strike + shift must be positive (otherwise the price does not depend on the vol)");
    const double intrinsic = std::max(optionSign(type) * (f - k), 0.0);
    const double tv = timeValue("blackImpliedVol", price, intrinsic);
    if (tv <= 0.0) return 0.0;

    // The out-of-the-money option carries the time value: a call when k ≥ f, a put otherwise. Its price rises
    // from 0 to its bound (f for the call, k for the put) as σ goes from 0 to ∞.
    const OptionType otm = k >= f ? OptionType::Call : OptionType::Put;
    const double bound = otm == OptionType::Call ? f : k;
    RL_REQUIRE(tv < bound * (1.0 - 1e-14), "blackImpliedVol: price is at or above the no-arbitrage bound");
    const double sqrtT = std::sqrt(expiry);
    auto otmPrice = [&](double sigma) { return blackPrice<double>(forward, strike, expiry, sigma, shift, otm); };

    // Bracket in σ: start from the Brenner–Subrahmanyam ATM estimate and expand geometrically.
    const double guess = std::clamp(tv * kSqrt2Pi / (std::sqrt(f * k) * sqrtT), 1e-8, 50.0);
    double lo = guess, hi = guess;
    for (int i = 0; i < 200 && otmPrice(lo) >= tv; ++i) lo *= 0.5;
    for (int i = 0; i < 200 && otmPrice(hi) <= tv; ++i) {
        lo = hi; // the previous upper end is now known to be below the root: keep the bracket a factor of 2 wide
        hi *= 2.0;
    }
    RL_REQUIRE(otmPrice(lo) < tv && otmPrice(hi) > tv, "blackImpliedVol: could not bracket the vol");
    if (!(lo < hi)) hi = lo * 2.0;

    auto fdf = [&](double sigma) -> std::pair<double, double> {
        const double p = otmPrice(sigma);
        const double v = blackVega<double>(forward, strike, expiry, sigma, shift);
        if (p <= kTiny) {
            const double xx = std::log(f / k);
            return {-800.0 - std::log(tv), xx * xx / (sigma * sigma * sigma * expiry)};
        }
        return {std::log(p / tv), v / p};
    };
    math::RootOptions ro;
    ro.xTol = options.relTol * hi;
    ro.maxIter = options.maxIter;
    try {
        // Start from the geometric middle of the bracket. For a far out-of-the-money option the ATM-style
        // guess is orders of magnitude too low and lands where the price underflows.
        return math::newtonSafeguarded(fdf, std::sqrt(lo * hi), lo, hi, ro).x;
    } catch (const NumericalError&) {
        // DESIGN: the fallback exists for the pathological flat-vega corners (σ√T ≫ 10); Brent on the price
        // itself cannot overshoot and needs no derivative.
        math::RootOptions rb;
        rb.xTol = options.relTol * hi;
        rb.maxIter = 500;
        return math::brent([&](double sigma) { return otmPrice(sigma) - tv; }, lo, hi, rb).x;
    }
}

} // namespace rl::vol
