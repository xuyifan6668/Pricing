// The model itself: its two functions H and ζ, the Hull–White mapping, and the property everything
// else rests on: deflated bonds are martingales, so the model reprices today's curve for any κ, ζ.
#include "rl/models/lgm.hpp"

#include "rl/core/error.hpp"
#include "rl/math/quadrature.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>

using rl::models::Lgm1F;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

TEST_CASE("H: (1 - exp(-kt))/k, with the k -> 0 limit H = t", "[models][lgm]") {
    CHECK_THAT(Lgm1F::H(0.03, 10.0), WithinRel((1.0 - std::exp(-0.3)) / 0.03, 1e-14));
    CHECK_THAT(Lgm1F::H(0.0, 7.0), WithinAbs(7.0, 1e-14));
    CHECK_THAT(Lgm1F::H(1e-10, 7.0), WithinRel(7.0, 1e-8));
    CHECK(Lgm1F::H(0.05, 0.0) == 0.0);
    // Negative mean reversion is legal (H grows faster than t): rates shocks amplify.
    CHECK(Lgm1F::H(-0.02, 10.0) > 10.0);
    // H is increasing and concave for k > 0: later cashflows are more exposed to the state, with diminishing effect.
    CHECK(Lgm1F::H(0.05, 20.0) - Lgm1F::H(0.05, 10.0) < Lgm1F::H(0.05, 10.0));
}

TEST_CASE("zeta: constant Hull-White sigma integrates to sigma^2 (exp(2kt) - 1)/(2k)", "[models][lgm]") {
    for (double kappa : {0.0, 0.03, 0.10, -0.02}) {
        const double sigma = 0.009;
        const Lgm1F m = Lgm1F::constantSigma(kappa, sigma, 30.0);
        for (double t : {0.25, 1.0, 5.0, 10.0, 30.0, 40.0}) {
            const double exact = kappa == 0.0 ? sigma * sigma * t : sigma * sigma * std::expm1(2.0 * kappa * t) / (2.0 * kappa);
            CHECK_THAT(m.zeta(t), WithinRel(exact, 1e-12));
        }
        CHECK(m.zeta(0.0) == 0.0);
        CHECK_THAT(m.hullWhiteSigmas().front(), WithinRel(sigma, 1e-12));
    }
}

TEST_CASE("zeta: reproduces its nodes, is increasing, and maps back to piecewise-constant sigma", "[models][lgm]") {
    const Lgm1F m(0.04, {1.0, 2.0, 5.0, 10.0}, {1.0e-4, 1.9e-4, 5.2e-4, 1.3e-3});
    for (std::size_t i = 0; i < 4; ++i) CHECK_THAT(m.zeta(m.zetaTimes()[i]), WithinRel(m.zetaValues()[i], 1e-13));
    double prev = 0.0;
    for (double t = 0.05; t < 12.0; t += 0.05) {
        CHECK(m.zeta(t) > prev);
        prev = m.zeta(t);
    }
    // Rebuilding zeta from the recovered sigmas gives the same nodes.
    const auto sig = m.hullWhiteSigmas();
    double z = 0.0, t0 = 0.0;
    for (std::size_t i = 0; i < 4; ++i) {
        const double t1 = m.zetaTimes()[i];
        z += sig[i] * sig[i] * (std::exp(0.08 * t1) - std::exp(0.08 * t0)) / 0.08;
        CHECK_THAT(z, WithinRel(m.zetaValues()[i], 1e-12));
        t0 = t1;
    }
}

TEST_CASE("Validation: zeta is an accumulated variance", "[models][lgm]") {
    CHECK_THROWS_AS(Lgm1F(0.03, {1.0, 2.0}, {2e-4, 1e-4}), rl::InputError); // decreasing
    CHECK_THROWS_AS(Lgm1F(0.03, {2.0, 1.0}, {1e-4, 2e-4}), rl::InputError); // times not increasing
    CHECK_THROWS_AS(Lgm1F(0.03, {1.0}, {0.0}), rl::InputError);
    CHECK_THROWS_AS(Lgm1F(0.03, {}, {}), rl::InputError);
    CHECK_NOTHROW(Lgm1F(0.03, {1.0, 2.0}, {1e-4, 1e-4})); // flat is allowed: zero vol on that interval
}

TEST_CASE("Martingale: E[P(t,T;x)/N(t,x)] = D(T) for every t < T, so the model reprices the curve", "[models][lgm][martingale]") {
    const rl::math::GaussHermite gh(48);
    for (double kappa : {0.0, 0.03, 0.15}) {
        const Lgm1F m = Lgm1F::constantSigma(kappa, 0.011, 30.0);
        for (double t : {0.5, 2.0, 10.0}) {
            for (double T : {t, t + 1.0, 30.0}) {
                const double dfT = std::exp(-0.037 * T);
                const double zeta = m.zeta(t);
                const double mean = gh.integrate(0.0, std::sqrt(zeta), [&](double x) { return m.deflatedBond(dfT, T, zeta, x); });
                CHECK_THAT(mean, WithinRel(dfT, 1e-12));
            }
        }
    }
}

TEST_CASE("Bond reconstruction: P(t,t) = 1, x = 0 at zero variance gives the forward bond, bonds fall as x rises", "[models][lgm]") {
    const Lgm1F m = Lgm1F::constantSigma(0.03, 0.01, 30.0);
    const double t = 3.0, T = 8.0, d0 = std::exp(-0.035 * t), d1 = std::exp(-0.035 * T);
    const double zeta = m.zeta(t);
    for (double x : {-0.05, 0.0, 0.04}) CHECK_THAT(m.bond(d0, t, d0, t, zeta, x), WithinAbs(1.0, 1e-14));
    CHECK_THAT(m.bond(d0, t, d1, T, 0.0, 0.0), WithinRel(d1 / d0, 1e-14));
    CHECK(m.bond(d0, t, d1, T, zeta, 0.02) < m.bond(d0, t, d1, T, zeta, 0.0));
    // The numeraire undoes the deflation.
    CHECK_THAT(m.deflatedBond(d1, T, zeta, 0.01) * m.numeraire(d0, t, zeta, 0.01), WithinRel(m.bond(d0, t, d1, T, zeta, 0.01), 1e-14));
}
