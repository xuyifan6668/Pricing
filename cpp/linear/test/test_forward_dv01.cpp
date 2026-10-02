#include "rl/linear/annuity.hpp"
#include "rl/linear/dv01.hpp"
#include "rl/linear/forward.hpp"
#include "rl/linear/leg_pv.hpp"
#include "rl/linear/par_rate.hpp"

#include "rl/aad/real.hpp"
#include "rl/core/error.hpp"
#include "rl/instruments/builder.hpp"
#include "test_support.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>

using namespace rl::linear;
using namespace rl::linear::test;
using rl::aad::Real;
using rl::core::Date;
using rl::core::DayCount;
using rl::core::Tenor;
using rl::instruments::OisSwapBuilder;
using rl::math::value;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

TEST_CASE("Forwards off a flat curve", "[linear][forward]") {
    const double r = 0.04;
    auto c = flatCurve<double>("c", r);
    CHECK_THAT(continuousForward(c, 1.0, 2.5), WithinRel(r, 1e-13));
    CHECK_THAT(simpleForward(c, 1.0, 2.0, 1.0), WithinRel(std::exp(r) - 1.0, 1e-13));
    const Date d1(2027, 3, 9), d2(2027, 6, 9);
    const double tau = rl::core::yearFraction(DayCount::Act360, d1, d2);
    const double dt = c.time(d2) - c.time(d1);
    CHECK_THAT(simpleForward(c, d1, d2, DayCount::Act360), WithinRel((std::exp(r * dt) - 1.0) / tau, 1e-13));
    CHECK_THROWS_AS(simpleForward(c, 2.0, 1.0, 1.0), rl::InputError);
    CHECK_THROWS_AS(continuousForward(c, 2.0, 2.0), rl::InputError);
}

TEST_CASE("Fixed-rate DV01 is sign * A * 1e-4; parallel DV01 has the right sign and size", "[linear][dv01]") {
    auto cs = flatSet<double>(0.035, 0.035);
    rl::market::FixingStore fx;
    auto payer = OisSwapBuilder::usdSofr(kAsOf).id("S").tenor(Tenor::years(10)).fixedRate(0.035).build();
    const double a = annuity(payer.fixedLeg(), cs);
    CHECK_THAT(fixedRateDv01(payer, cs), WithinRel(-a * 1e-4, 1e-14)); // payer: PV falls as the fixed rate rises
    const double pd = parallelDv01(payer, cs, fx);
    CHECK(pd > 0.0);                          // payer gains when rates rise
    CHECK_THAT(pd, WithinRel(a * 1e-4, 0.05)); // near par: ≈ annuity per bp
    // Second-order: the +1bp and −1bp one-sided bumps differ by convexity × 1bp ≈ 0.1% of DV01 on a 10y.
    const double pdDown = -parallelDv01(payer, cs, fx, -1.0);
    CHECK_THAT(pd, WithinRel(pdDown, 1e-2));
    CHECK(pd - pdDown < 0.0); // the payer is short the fixed-rate bond: negative convexity, gains less on the way up
}

TEST_CASE("parallelBump shifts every curve's zero rate by exactly the bump", "[linear][dv01]") {
    auto cs = flatSet<double>(0.035, 0.037);
    auto up = parallelBump(cs, 0.0010);
    for (const auto& name : cs.names()) {
        const auto& c0 = cs.curve(name);
        const auto& c1 = up.curve(name);
        for (double t : {0.3, 1.0, 4.0, 12.0, 25.0}) CHECK_THAT(-std::log(c1.df(t)) / t, WithinAbs(-std::log(c0.df(t)) / t + 0.0010, 1e-13));
    }
    CHECK(up.projection("USD-SOFR").df(1.0) == up.discount().df(1.0)); // roles survive the bump
}

TEST_CASE("DV01 under Real matches double", "[linear][dv01][aad]") {
    auto cd = flatSet<double>(0.035, 0.035);
    auto cr = flatSet<Real>(0.035, 0.035);
    rl::market::FixingStore fx;
    auto s = OisSwapBuilder::usdSofr(kAsOf).id("S").tenor(Tenor::years(7)).fixedRate(0.03).build();
    CHECK(fixedRateDv01(s, cd) == value(fixedRateDv01(s, cr)));
    CHECK(parallelDv01(s, cd, fx) == value(parallelDv01(s, cr, fx)));
}
