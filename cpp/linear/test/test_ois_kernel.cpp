// The overnight compounding kernel: projected, realised, mixed, lookback, lockout, averaged, missing fixings.
#include "rl/linear/leg_pv.hpp"

#include "rl/aad/real.hpp"
#include "rl/core/error.hpp"
#include "rl/instruments/builder.hpp"
#include "test_support.hpp"

#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>

using namespace rl::linear;
using namespace rl::linear::test;
using rl::aad::Real;
using rl::core::Calendar;
using rl::core::Date;
using rl::core::Index;
using rl::core::Tenor;
using rl::instruments::Side;
using rl::math::value;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

TEMPLATE_TEST_CASE("Fully projected coupon: compounded rate is the DF ratio over the observation period", "[linear][ois]", double, Real) {
    using T = TestType;
    const double r = 0.0385;
    auto cs = flatSet<T>(r, r);
    rl::market::FixingStore fx; // no fixings needed
    auto c = oisCoupon(Date(2027, 3, 9), Date(2027, 6, 9), Date(2027, 6, 11), 1e6);
    const auto& proj = cs.projection("USD-SOFR");
    const double expected = value((proj.df(Date(2027, 3, 9)) / proj.df(Date(2027, 6, 9)) - 1.0)) / c.tau;
    CHECK_THAT(value(overnightRate(c, cs, fx)), WithinRel(expected, 1e-14));
    CHECK_THAT(value(pv(c, cs, fx)), WithinRel(1e6 * expected * c.tau * value(cs.discount().df(Date(2027, 6, 11))), 1e-14));
}

TEST_CASE("A spread-free OIS leg with payment on accrual end telescopes to P(start) - P(end)", "[linear][ois]") {
    // Any curve shape; the compounded coupon on the projection curve, discounted on the same curve, telescopes.
    std::vector<double> ts{0.5, 1.0, 2.0, 5.0, 10.0, 30.0};
    std::vector<double> ys{0.019, 0.035, 0.066, 0.17, 0.36, 1.1};
    rl::market::YieldCurve<double> c(kAsOf, "SOFR", ts, ys, rl::market::CurveInterpolation::LinearLogDf); // non-flat; MonotoneConvex works too
    rl::market::CurveSet<double> cs({c}, "SOFR", {{"USD-SOFR", "SOFR"}});
    rl::market::FixingStore fx;
    rl::core::ScheduleSpec spec;
    spec.effective = Date(2026, 9, 9);
    spec.termination = Date(2031, 9, 9);
    spec.frequency = Tenor::months(6);
    spec.calendar = Calendar::get("USD-SIFMA");
    spec.paymentLag = 0;
    auto sched = rl::core::Schedule::generate(spec);
    auto leg = rl::instruments::makeOvernightLeg(sched, Index::sofr(), 0.0, 1.0, Side::Receive, rl::core::Currency::USD);
    const double expected = c.df(sched.front().start) - c.df(sched.back().end);
    CHECK_THAT(pv(leg, cs, fx), WithinAbs(expected, 1e-14));
}

TEST_CASE("Partially elapsed coupon: realised fixings compound into the projected remainder", "[linear][ois]") {
    const double r = 0.0385;
    auto cs = flatSet<double>(r, r);
    const Index sofr = Index::sofr();
    const Calendar cal = Calendar::get(sofr.fixingCalendar);
    const double rf = 0.0400; // realised fixings above the curve
    auto fx = constantFixings(sofr, Date(2026, 8, 1), Date(2026, 9, 3), rf);
    // Coupon started 2026-08-10, ends 2026-11-10; as-of 2026-09-04 (Friday) the fixings for 10 Aug ... 3 Sep are known.
    auto c = oisCoupon(Date(2026, 8, 10), Date(2026, 11, 10), Date(2026, 11, 12), 1e6);
    double growth = 1.0;
    for (Date d(2026, 8, 10); d < kAsOf; d = cal.advance(d, 1)) growth *= 1.0 + rf * rl::core::yearFraction(sofr.dayCount, d, cal.advance(d, 1));
    const auto& proj = cs.projection("USD-SOFR");
    growth *= proj.df(kAsOf) / proj.df(Date(2026, 11, 10));
    const double expected = (growth - 1.0) / c.tau;
    CHECK_THAT(overnightRate(c, cs, fx), WithinRel(expected, 1e-13));
    // The realised part is above the curve, so the rate sits between the two.
    CHECK(overnightRate(c, cs, fx) > r);
    CHECK(overnightRate(c, cs, fx) < rf + 1e-4);
    // Fully realised coupon: no curve involved at all.
    auto past = oisCoupon(Date(2026, 8, 3), Date(2026, 9, 1), Date(2026, 9, 3), 1e6);
    double g2 = 1.0;
    for (Date d(2026, 8, 3); d < Date(2026, 9, 1); d = cal.advance(d, 1)) g2 *= 1.0 + rf * rl::core::yearFraction(sofr.dayCount, d, cal.advance(d, 1));
    CHECK_THAT(overnightRate(past, cs, fx), WithinRel((g2 - 1.0) / past.tau, 1e-13));
    CHECK(pv(past, cs, fx) == 0.0); // paid 3 Sep, before as-of
}

TEST_CASE("Missing fixings throw with the index, date and coupon in the message", "[linear][ois]") {
    auto cs = flatSet<double>(0.0385, 0.0385);
    rl::market::FixingStore fx;
    fx.add("USD-SOFR", Date(2026, 9, 1), 0.04);
    auto c = oisCoupon(Date(2026, 8, 31), Date(2026, 11, 30), Date(2026, 12, 2), 1e6);
    try {
        overnightRate(c, cs, fx);
        FAIL("expected InputError");
    } catch (const rl::InputError& e) {
        const std::string w = e.what();
        CHECK(w.find("USD-SOFR") != std::string::npos);
        CHECK(w.find("2026-08-31") != std::string::npos);
        CHECK(w.find("2026-11-30") != std::string::npos);
    }
}

TEST_CASE("Lookback shifts the observation window; lockout freezes the last fixings", "[linear][ois]") {
    const double r = 0.0385;
    auto cs = flatSet<double>(r, r);
    const Index sofr = Index::sofr();
    const Calendar cal = Calendar::get(sofr.fixingCalendar);
    // Distinct fixings per day so that which day is used is observable.
    rl::market::FixingStore fx;
    int k = 0;
    for (Date d(2026, 7, 1); d <= Date(2026, 9, 3); ++d)
        if (cal.isBusinessDay(d)) fx.add(sofr.name, d, 0.03 + 0.0001 * (k++));

    SECTION("lookback 5 on a fully realised coupon equals the plain rate on the shifted dates") {
        auto c = oisCoupon(Date(2026, 8, 3), Date(2026, 9, 1), Date(2026, 9, 3), 1e6, sofr, 0.0, 5, 0);
        auto shifted = oisCoupon(cal.advance(Date(2026, 8, 3), -5), cal.advance(Date(2026, 9, 1), -5), Date(2026, 9, 3), 1e6);
        // Same compounding factor; different tau (observation vs accrual) — compare the growth factors.
        const double gL = 1.0 + overnightRate(c, cs, fx) * rl::core::yearFraction(sofr.dayCount, shifted.accrualStart, shifted.accrualEnd);
        const double gS = 1.0 + overnightRate(shifted, cs, fx) * shifted.tau;
        CHECK_THAT(gL, WithinRel(gS, 1e-14));
    }
    SECTION("lockout 2 uses the fixing of the second-last business day for the last two days") {
        const Date s(2026, 8, 3), e(2026, 9, 1);
        auto locked = oisCoupon(s, e, Date(2026, 9, 3), 1e6, sofr, 0.0, 0, 2);
        const Date frozen = cal.advance(e, -2);
        double g = 1.0;
        for (Date d = s; d < e; d = cal.advance(d, 1)) {
            const Date fd = d >= frozen ? frozen : d;
            g *= 1.0 + fx.require(sofr.name, fd) * rl::core::yearFraction(sofr.dayCount, d, cal.advance(d, 1));
        }
        CHECK_THAT(overnightRate(locked, cs, fx), WithinRel((g - 1.0) / locked.tau, 1e-13));
        auto plain = oisCoupon(s, e, Date(2026, 9, 3), 1e6);
        CHECK(overnightRate(locked, cs, fx) != overnightRate(plain, cs, fx));
    }
}

TEST_CASE("Averaged index (Fed Funds) uses the arithmetic average; compounded SOFR differs by the compounding basis", "[linear][ois]") {
    const double r = 0.0385;
    auto cs = flatSet<double>(r, r);
    rl::market::FixingStore fx;
    auto ff = oisCoupon(Date(2027, 3, 9), Date(2027, 6, 9), Date(2027, 6, 11), 1e6, Index::fedFunds());
    auto so = oisCoupon(Date(2027, 3, 9), Date(2027, 6, 9), Date(2027, 6, 11), 1e6, Index::sofr());
    // On a flat continuously-compounded curve the average rate is r·(t2−t1)/τ (curve time over ACT/360 time).
    const auto& proj = cs.projection("USD-FEDFUNDS");
    const double dt = proj.time(Date(2027, 6, 9)) - proj.time(Date(2027, 3, 9));
    CHECK_THAT(overnightRate(ff, cs, fx), WithinRel(r * dt / ff.tau, 1e-13));
    // Compounded − averaged ≈ r² τ / 2 ≈ 1.9 bp on a quarter at 3.85%.
    const double basis = overnightRate(so, cs, fx) - overnightRate(ff, cs, fx);
    CHECK(basis > 1.0e-4);
    CHECK(basis < 3.0e-4);
    CHECK_THAT(basis, WithinRel(0.5 * r * r * dt * dt / ff.tau, 2e-2));
}

TEST_CASE("Kernel: Real equals double bit for bit, including realised fixings", "[linear][ois][aad]") {
    const double r = 0.0385;
    auto cd = flatSet<double>(r, r);
    std::vector<double> ts{0.5, 1.0, 2.0, 5.0, 10.0, 30.0};
    std::vector<Real> ys;
    for (double t : ts) ys.push_back(Real::input(r * t));
    rl::market::CurveSet<Real> cr({rl::market::YieldCurve<Real>(kAsOf, "USD-SOFR-OIS", ts, ys, rl::market::CurveInterpolation::LinearLogDf)}, "USD-SOFR-OIS",
                                  {{"USD-SOFR", "USD-SOFR-OIS"}});
    auto fx = constantFixings(Index::sofr(), Date(2026, 8, 1), Date(2026, 9, 3), 0.04);
    auto c = oisCoupon(Date(2026, 8, 10), Date(2026, 11, 10), Date(2026, 11, 12), 1e6);
    CHECK(overnightRate(c, cd, fx) == overnightRate(c, cr, fx).value());
    CHECK(overnightRate(c, cr, fx).isActive());
    rl::aad::Tape::current().reset();
}
