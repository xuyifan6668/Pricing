// Arbitrage diagnostic: Bachelier density sanity (mass, mean, peak), a benign node is clean, the stressed node is
// flagged at a grid-independent strike, unusable vols and call-spread violations are reported.
#include "rl/vol/arbitrage.hpp"

#include "rl/core/error.hpp"
#include "rl/math/normal.hpp"
#include "test_support.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <limits>

using namespace rl::vol;
using rl::market::SabrParams;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

TEST_CASE("strikeGrid and diagnosticStrikeRange", "[vol][arbitrage]") {
    const auto g = strikeGrid(0.01, 0.05, 5);
    REQUIRE(g.size() == 5);
    CHECK(g.front() == 0.01);
    CHECK(g.back() == 0.05);
    CHECK_THAT(g[2], WithinAbs(0.03, 1e-17));
    CHECK_THROWS_AS(strikeGrid(0.05, 0.01, 5), rl::InputError);
    CHECK_THROWS_AS(strikeGrid(0.01, 0.05, 1), rl::InputError);
    const auto p = test::kStressed.params();
    const auto [lo, hi] = diagnosticStrikeRange(p, test::kStressed.forward, test::kStressed.expiry, test::kStressed.atm);
    CHECK(lo == 1e-4); // F − 4σ√T < 0 on a 10y expiry: floored for β > 0
    CHECK_THAT(hi, WithinRel(test::kStressed.forward + 4.0 * test::kStressed.atm * std::sqrt(10.0), 1e-15));
}

TEST_CASE("Constant normal vol (beta = 0, nu = 0): the density is the normal density", "[vol][arbitrage]") {
    SabrParams<double> p;
    p.alpha = 0.01;
    p.beta = 0.0;
    p.rho = 0.0;
    p.nu = 0.0;
    const double f = 0.03, tt = 4.0, sd = 0.01 * 2.0;
    const auto d = riskNeutralDensity(p, f, tt, strikeGrid(f - 8.0 * sd, f + 8.0 * sd, 4001));
    const auto c = checkDensity(d, f);
    CHECK(c.arbitrageFree());
    CHECK_THAT(c.integratedMass, WithinAbs(1.0, 2e-6));
    CHECK_THAT(c.impliedMean, WithinAbs(f, 1e-7));
    const auto mid = d[2000];
    CHECK(mid.strike == f);
    CHECK_THAT(mid.density, WithinRel(rl::math::normalPdf(0.0) / sd, 1e-5));
    CHECK_THAT(mid.dCalldK, WithinAbs(-0.5, 1e-9));
    CHECK(c.minDensity > -1e-6); // the second difference's rounding noise in the far tails (arbitrage.hpp), not arbitrage
}

TEST_CASE("A benign 1Yx10Y node is arbitrage-free across four standard deviations", "[vol][arbitrage]") {
    const auto& n = test::kBenign;
    const auto p = n.params();
    const auto [lo, hi] = diagnosticStrikeRange(p, n.forward, n.expiry, n.atm, 4.0);
    const auto c = checkSabrSmile(p, n.forward, n.expiry, lo, hi, 801);
    CHECK(c.arbitrageFree());
    CHECK(c.minDensity > 0.0);
    CHECK_FALSE(c.firstNegativeDensityBelowForward.has_value());
}

TEST_CASE("The stressed 10Yx30Y node goes negative ~170 bp below the forward, independent of the grid", "[vol][arbitrage]") {
    const auto& n = test::kStressed;
    const auto p = n.params();
    const auto [lo, hi] = diagnosticStrikeRange(p, n.forward, n.expiry, n.atm, 4.0);
    const auto coarse = checkSabrSmile(p, n.forward, n.expiry, lo, hi, 401);
    const auto fine = checkSabrSmile(p, n.forward, n.expiry, lo, hi, 4001);
    REQUIRE(coarse.firstNegativeDensityBelowForward.has_value());
    REQUIRE(fine.firstNegativeDensityBelowForward.has_value());
    const double step = (hi - lo) / 400.0;
    CHECK(std::fabs(*coarse.firstNegativeDensityBelowForward - *fine.firstNegativeDensityBelowForward) <= step);
    CHECK(*fine.firstNegativeDensityBelowForward > 0.019);
    CHECK(*fine.firstNegativeDensityBelowForward < 0.023);
    CHECK_FALSE(fine.arbitrageFree());
    CHECK(fine.minDensity < 0.0);
    // A fit reprices the quoted strikes (F ± 100 bp) perfectly; the density is still fine there. The problem is below.
    CHECK(*fine.firstNegativeDensityBelowForward < n.forward - 0.01);
}

TEST_CASE("Unusable vols and call-spread violations are reported", "[vol][arbitrage]") {
    // β = 0.5, large α, 30y: the −β(2−β)α²/(24 √(FK)) term drives the vol below zero near K = 0.
    SabrParams<double> p;
    p.alpha = 0.08;
    p.beta = 0.5;
    p.rho = 0.0;
    p.nu = 0.0;
    const auto d = riskNeutralDensity(p, 0.03, 30.0, {1e-4, 0.02, 0.03, 0.04}, 1e-5);
    CHECK(d[0].normalVol < 0.0);
    CHECK(std::isnan(d[0].density));
    const auto c = checkDensity(d, 0.03);
    CHECK(c.nonPositiveVol);
    // The scan stops at the bad strike nearest the forward: on this smile the density is already negative at 2%,
    // far above where the vol itself goes negative — which is exactly why the report scans from the forward.
    CHECK(c.firstNegativeDensityBelowForward.has_value());
    CHECK_FALSE(c.arbitrageFree());
    CHECK_THROWS_AS(riskNeutralDensity(p, 0.03, 30.0, {1e-6}, 1e-5), rl::InputError); // K − h + s ≤ 0

    const double nan = std::numeric_limits<double>::quiet_NaN();
    std::vector<DensityPoint> manual{{0.02, 0.01, 0.011, -0.9, 10.0}, {0.03, 0.01, 0.004, 0.1, 10.0}, {0.04, 0.01, 0.001, -0.2, 10.0}};
    CHECK(checkDensity(manual, 0.03).callSpreadViolated);
    std::vector<DensityPoint> unsorted{{0.03, 0.01, nan, nan, nan}, {0.02, 0.01, nan, nan, nan}};
    CHECK_THROWS_AS(checkDensity(unsorted, 0.03), rl::InputError);

    // Scan order on a hand-built grid: nearest bad strike to the forward on each side; an unusable vol counts as bad.
    std::vector<DensityPoint> scan{{0.01, -0.001, nan, nan, nan}, {0.02, 0.01, 0.02, -0.5, 5.0}, {0.03, 0.01, 0.01, -0.5, 10.0},
                                   {0.04, 0.01, 0.005, -0.3, -1.0}, {0.05, 0.01, 0.002, -0.1, 3.0}};
    const auto s = checkDensity(scan, 0.03);
    REQUIRE(s.firstNegativeDensityBelowForward.has_value());
    CHECK(*s.firstNegativeDensityBelowForward == 0.01);
    REQUIRE(s.firstNegativeDensityAboveForward.has_value());
    CHECK(*s.firstNegativeDensityAboveForward == 0.04);
    CHECK(s.nonPositiveVol);
    CHECK(s.minDensity == -1.0);
    CHECK(s.minDensityStrike == 0.04);
    scan[1].density = -1.0;
    CHECK(*checkDensity(scan, 0.03).firstNegativeDensityBelowForward == 0.02);
}
