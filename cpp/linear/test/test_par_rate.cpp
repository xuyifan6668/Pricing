#include "rl/linear/annuity.hpp"
#include "rl/linear/leg_pv.hpp"
#include "rl/linear/par_rate.hpp"

#include "rl/core/error.hpp"
#include "rl/instruments/builder.hpp"
#include "test_support.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>

using namespace rl::linear;
using namespace rl::linear::test;
using rl::core::Date;
using rl::core::Tenor;
using rl::instruments::BasisSwapBuilder;
using rl::instruments::OisSwapBuilder;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

TEST_CASE("Annuity is the sum of tau * N * P(pay) over unsettled coupons", "[linear][annuity]") {
    const double r = 0.04;
    auto cs = flatSet<double>(r, r);
    auto s = OisSwapBuilder::usdSofr(kAsOf).id("S").tenor(Tenor::years(3)).fixedRate(0.03).notional(1.0).build();
    double a = 0.0;
    for (const auto& cf : s.fixedLeg().cashflows()) {
        const auto& c = std::get<rl::instruments::FixedCoupon>(cf);
        a += c.tau * cs.discount().df(c.payment);
    }
    CHECK_THAT(annuity(s.fixedLeg(), cs), WithinRel(a, 1e-14));
    CHECK(annuity(s.fixedLeg(), cs) == annuity(s.floatLeg(), cs)); // same schedule, same τ convention
    // A 3y annual annuity at 4% is a bit under 3 × e^{-0.08}.
    CHECK(a > 2.6);
    CHECK(a < 2.9);
}

TEST_CASE("Par rate zeroes the swap and satisfies A * K = -sign * PV(float)", "[linear][par]") {
    auto cs = flatSet<double>(0.035, 0.035);
    rl::market::FixingStore fx;
    for (int years : {1, 2, 5, 10, 30}) {
        auto b = OisSwapBuilder::usdSofr(kAsOf).id("S").tenor(Tenor::years(years)).fixedRate(0.03);
        auto s = b.build();
        const double k = parRate(s, cs, fx);
        auto par = b.fixedRate(k).build();
        CHECK_THAT(pv(par, cs, fx), WithinAbs(0.0, 1e-8)); // 1e-14 of a 1e6 notional
        CHECK_THAT(annuity(s.fixedLeg(), cs) * k, WithinRel(-s.fixedLeg().sign() * pv(s.floatLeg(), cs, fx), 1e-13));
        // On a flat 3.5% continuously-compounded curve the annual-pay ACT/360 par rate is a little above 3.5%.
        CHECK(k > 0.035);
        CHECK(k < 0.037);
        // Receiver and payer agree on the par rate.
        CHECK_THAT(parRate(b.payFixed(false).build(), cs, fx), WithinRel(k, 1e-13));
    }
}

TEST_CASE("Par spread zeroes the basis swap; identical curves give the compounding basis", "[linear][par][basis]") {
    const double r = 0.0385;
    auto cs = flatSet<double>(r, r);
    rl::market::FixingStore fx;
    auto b = BasisSwapBuilder::usdFedFundsSofr(kAsOf).id("B").tenor(Tenor::years(5)).spread(-0.0002);
    auto bs = b.build();
    const double s = parSpread(bs, cs, fx, 0);
    // FF averaged vs SOFR compounded on the same curve: FF leg must receive the compounding basis, ~+1.9 bp.
    CHECK(s > 1.0e-4);
    CHECK(s < 3.0e-4);
    auto par = b.spread(s).build();
    CHECK_THAT(pv(par, cs, fx), WithinAbs(0.0, 1e-8));
    // Solving for the spread on the flat leg instead: same schedule, so A_1 = A_0 and s_1 = s_cur − s_0 exactly.
    CHECK_THAT(parSpread(bs, cs, fx, 1), WithinRel(-0.0002 - s, 1e-10));
    CHECK_THROWS_AS(parSpread(bs, cs, fx, 2), rl::InputError);
    // A lower FF curve needs a higher FF spread.
    auto cs2 = flatSet<double>(r, r - 0.0010);
    CHECK(parSpread(bs, cs2, fx, 0) > s + 0.0009);
}
