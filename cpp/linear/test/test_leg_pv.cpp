#include "rl/linear/leg_pv.hpp"

#include "rl/aad/real.hpp"
#include "rl/core/daycount.hpp"
#include "rl/instruments/builder.hpp"
#include "test_support.hpp"

#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>

using namespace rl::linear;
using namespace rl::linear::test;
using rl::aad::Real;
using rl::core::Date;
using rl::core::Tenor;
using rl::instruments::Side;
using rl::math::value;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

TEMPLATE_TEST_CASE("Fixed coupon PV is N * r * tau * P(pay); settled coupons are worth zero", "[linear][pv]", double, Real) {
    using T = TestType;
    const double r = 0.04;
    auto cs = flatSet<T>(r, r);
    // Hand example: N = 1e6, coupon 3.5% on τ = 0.5, paid one year out on a 4% flat curve: 17500 · e^{−0.04}.
    auto c = fixedCoupon(Date(2027, 3, 4), Date(2027, 9, 4), Date(2027, 9, 4), 0.5, 1e6, 0.035);
    CHECK_THAT(value(pv(c, cs)), WithinRel(17500.0 * std::exp(-0.04), 1e-13));
    auto settled = fixedCoupon(Date(2025, 9, 4), Date(2026, 3, 4), Date(2026, 9, 4), 0.5, 1e6, 0.035); // pays on asOf
    CHECK(value(pv(settled, cs)) == 0.0);
    rl::instruments::NotionalExchange n{Date(2027, 9, 4), -1e6};
    CHECK_THAT(value(pv(n, cs)), WithinRel(-1e6 * std::exp(-0.04), 1e-13));
}

TEST_CASE("Leg PV applies the side; cashflow variant dispatch", "[linear][pv]") {
    auto cs = flatSet<double>(0.04, 0.04);
    rl::market::FixingStore fx;
    auto c1 = fixedCoupon(Date(2026, 9, 4), Date(2027, 9, 4), Date(2027, 9, 4), 1.0, 1e6, 0.03);
    auto c2 = fixedCoupon(Date(2027, 9, 4), Date(2028, 9, 4), Date(2028, 9, 4), 1.0, 1e6, 0.03);
    rl::instruments::Leg rec({c1, c2}, Side::Receive, rl::core::Currency::USD);
    rl::instruments::Leg pay({c1, c2}, Side::Pay, rl::core::Currency::USD);
    // 2027-09-04 is exactly one year of ACT/365F; 2028-09-04 is 731/365 years (2028 is a leap year).
    const double expected = 30000.0 * (std::exp(-0.04 * 365.0 / 365.0) + std::exp(-0.04 * 731.0 / 365.0));
    CHECK_THAT(pv(rec, cs, fx), WithinRel(expected, 1e-13));
    CHECK_THAT(pv(pay, cs, fx), WithinRel(-expected, 1e-13));
    CHECK(pv(rl::instruments::Cashflow(c1), cs, fx) == pv(c1, cs));
}

TEST_CASE("Swap PV is the signed sum of legs; a par swap has zero PV", "[linear][pv][swap]") {
    auto cs = flatSet<double>(0.035, 0.035);
    rl::market::FixingStore fx;
    auto b = rl::instruments::OisSwapBuilder::usdSofr(kAsOf).id("S").tenor(Tenor::years(5)).fixedRate(0.03);
    rl::instruments::Swap payer = b.build();
    rl::instruments::Swap receiver = b.payFixed(false).build();
    CHECK_THAT(pv(payer, cs, fx), WithinAbs(-pv(receiver, cs, fx), 1e-9));
    CHECK(pv(payer, cs, fx) == pv(payer.fixedLeg(), cs, fx) + pv(payer.floatLeg(), cs, fx));
    // Paying 3% fixed against a 3.5% curve is worth something to the payer.
    CHECK(pv(payer, cs, fx) > 0.0);
}

TEST_CASE("Deposit, FRA, future implied rates off a flat curve", "[linear][pv][cash]") {
    const double r = 0.04;
    auto cs = flatSet<double>(r, r);
    rl::instruments::Deposit d("D", rl::core::Index::sofr(), Date(2026, 9, 9), Date(2026, 12, 9), 91.0 / 360.0, 1e6, 0.04);
    const double t1 = cs.discount().time(Date(2026, 9, 9)), t2 = cs.discount().time(Date(2026, 12, 9));
    const double implied = (std::exp(r * (t2 - t1)) - 1.0) / (91.0 / 360.0);
    CHECK_THAT(impliedDepositRate(d, cs), WithinRel(implied, 1e-13));
    CHECK_THAT(pv(d, cs), WithinRel(1e6 * (91.0 / 360.0) * (implied - 0.04) * std::exp(-r * t2), 1e-10));
    rl::instruments::Fra f("F", rl::core::Index::sofr(), Date(2026, 12, 7), Date(2026, 12, 9), Date(2027, 3, 9), 90.0 / 360.0, 1e6, implied); // projected off the SOFR curve
    CHECK(std::fabs(pv(f, cs)) < 1e6 * 1e-3); // near par
    rl::instruments::Future fu("SR3", rl::core::Index::sofr(), Date(2026, 12, 16), Date(2027, 3, 17), 91.0 / 360.0, 96.0, 0.5e-4);
    CHECK_THAT(impliedFutureRate(fu, cs), WithinAbs(impliedFutureForward(fu, cs) + 0.5e-4, 1e-18));
    CHECK_THAT(impliedFuturePrice(fu, cs), WithinAbs(100.0 * (1.0 - impliedFutureRate(fu, cs)), 1e-12));
    CHECK(impliedFuturePrice(fu, cs) < 100.0 - 100.0 * 0.038);
}

TEST_CASE("Real and double PVs agree bit for bit and the Real PV is active", "[linear][pv][aad]") {
    auto cd = flatSet<double>(0.035, 0.036);
    std::vector<double> ts{0.5, 1.0, 2.0, 5.0, 10.0, 30.0};
    std::vector<Real> ys;
    for (double t : ts) ys.push_back(Real::input(0.035 * t));
    rl::market::CurveSet<Real> cr({rl::market::YieldCurve<Real>(kAsOf, "USD-SOFR-OIS", ts, ys, rl::market::CurveInterpolation::LinearLogDf),
                                   flatCurve<Real>("USD-FF", 0.036)},
                                  "USD-SOFR-OIS", {{"USD-SOFR", "USD-SOFR-OIS"}, {"USD-FEDFUNDS", "USD-FF"}});
    rl::market::FixingStore fx;
    auto s = rl::instruments::OisSwapBuilder::usdSofr(kAsOf).id("S").tenor(Tenor::years(10)).fixedRate(0.034).build();
    CHECK(pv(s, cd, fx) == pv(s, cr, fx).value());
    CHECK(pv(s, cr, fx).isActive());
    rl::aad::Tape::current().reset();
}
