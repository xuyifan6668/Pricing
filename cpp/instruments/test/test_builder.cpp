// Trade builders.
#include "rl/instruments/builder.hpp"

#include "rl/core/error.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <variant>

using namespace rl::instruments;
using rl::core::Date;
using rl::core::Tenor;
using Catch::Matchers::WithinAbs;

namespace {
const Date kTrade(2026, 9, 4); // Friday; Labor Day Monday 7 Sep → spot Wednesday 9 Sep

template <class C>
const C& as(const Cashflow& cf) {
    const C* p = std::get_if<C>(&cf);
    REQUIRE(p != nullptr);
    return *p;
}
} // namespace

TEST_CASE("OisSwapBuilder: 10y USD SOFR from an explicit effective date matches the hand-written schedule", "[instruments][builder]") {
    // Same dates as cpp/core/test/test_schedule.cpp ("10y annual SOFR schedule from 2026-09-08").
    Swap s = OisSwapBuilder::usdSofr(kTrade).id("SOFR-10Y").effective(Date(2026, 9, 8)).tenor(Tenor::years(10)).fixedRate(0.0345).build();
    CHECK(s.id() == "SOFR-10Y");
    REQUIRE(s.fixedLeg().cashflows().size() == 10);
    REQUIRE(s.floatLeg().cashflows().size() == 10);
    CHECK(s.fixedLeg().sign() == -1.0); // payer by default
    CHECK(s.floatLeg().sign() == 1.0);
    const auto& c3 = as<FixedCoupon>(s.fixedLeg().cashflows()[3]);
    CHECK(c3.accrualStart == Date(2029, 9, 10));
    CHECK(c3.accrualEnd == Date(2030, 9, 9));
    CHECK(c3.payment == Date(2030, 9, 11));
    CHECK_THAT(c3.tau, WithinAbs(364.0 / 360.0, 1e-15)); // ACT/360, 2029-09-10 → 2030-09-09 is 364 days
    CHECK(c3.rate == 0.0345);
    CHECK(c3.notional == 1.0e6);
    const auto& f0 = as<OvernightCoupon>(s.floatLeg().cashflows()[0]);
    CHECK(f0.index.name == "USD-SOFR");
    CHECK(f0.accrualStart == Date(2026, 9, 8));
    CHECK(f0.payment == Date(2027, 9, 10));
    CHECK(f0.spread == 0.0);
    CHECK(s.maturity() == Date(2036, 9, 10));
}

TEST_CASE("OisSwapBuilder: spot rule, short tenors are single-period, options", "[instruments][builder]") {
    OisSwapBuilder b = OisSwapBuilder::usdSofr(kTrade);
    CHECK(b.effectiveDate() == Date(2026, 9, 9));
    Swap w = b.id("1W").tenor(Tenor::weeks(1)).fixedRate(0.0385).build();
    REQUIRE(w.fixedLeg().cashflows().size() == 1);
    const auto& c = as<FixedCoupon>(w.fixedLeg().cashflows()[0]);
    CHECK(c.accrualStart == Date(2026, 9, 9));
    CHECK(c.accrualEnd == Date(2026, 9, 16));
    CHECK(c.payment == Date(2026, 9, 18));
    CHECK_THAT(c.tau, WithinAbs(7.0 / 360.0, 1e-16));
    Swap y1 = OisSwapBuilder::usdSofr(kTrade).id("1Y").tenor(Tenor::years(1)).fixedRate(0.035).build();
    CHECK(y1.fixedLeg().cashflows().size() == 1); // ≤ 1Y: one period
    Swap m18 = OisSwapBuilder::usdSofr(kTrade).id("18M").tenor(Tenor::months(18)).fixedRate(0.035).build();
    CHECK(m18.fixedLeg().cashflows().size() == 2); // short front stub + one annual period
    CHECK(as<FixedCoupon>(m18.fixedLeg().cashflows()[0]).accrualEnd == Date(2027, 3, 9));
    Swap rec = OisSwapBuilder::usdSofr(kTrade).id("R").tenor(Tenor::years(2)).fixedRate(0.03).payFixed(false).notional(5e6).floatSpread(0.0001).lookback(5).build();
    CHECK(rec.fixedLeg().sign() == 1.0);
    CHECK(rec.floatLeg().sign() == -1.0);
    CHECK(as<FixedCoupon>(rec.fixedLeg().cashflows()[0]).notional == 5e6);
    CHECK(as<OvernightCoupon>(rec.floatLeg().cashflows()[0]).spread == 0.0001);
    CHECK(as<OvernightCoupon>(rec.floatLeg().cashflows()[0]).lookbackDays == 5);
    CHECK_THROWS_AS(OisSwapBuilder::usdSofr(kTrade).tenor(Tenor::years(1)).fixedRate(0.03).build(), rl::InputError); // no id
    CHECK_THROWS_AS(OisSwapBuilder::usdSofr(kTrade).id("x").fixedRate(0.03).build(), rl::InputError);               // no tenor
    CHECK_THROWS_AS(OisSwapBuilder::usdSofr(kTrade).id("x").tenor(Tenor::years(1)).build(), rl::InputError);         // no rate
}

TEST_CASE("BasisSwapBuilder: FF+spread vs SOFR, quarterly, spread leg first", "[instruments][builder]") {
    BasisSwap b = BasisSwapBuilder::usdFedFundsSofr(kTrade).id("FF-SOFR-2Y").tenor(Tenor::years(2)).spread(-0.00025).build();
    REQUIRE(b.legs()[0].cashflows().size() == 8);
    REQUIRE(b.legs()[1].cashflows().size() == 8);
    const auto& ff = as<OvernightCoupon>(b.legs()[0].cashflows()[0]);
    const auto& so = as<OvernightCoupon>(b.legs()[1].cashflows()[0]);
    CHECK(ff.index.name == "USD-FEDFUNDS");
    CHECK(ff.spread == -0.00025);
    CHECK(so.index.name == "USD-SOFR");
    CHECK(so.spread == 0.0);
    CHECK(ff.accrualStart == so.accrualStart);
    CHECK(ff.payment == so.payment);
    CHECK(ff.accrualStart == Date(2026, 9, 9));
    CHECK(ff.accrualEnd == Date(2026, 12, 9));
    CHECK(b.legs()[0].sign() == -1.0);
    CHECK(b.legs()[1].sign() == 1.0);
    BasisSwap m3 = BasisSwapBuilder::usdFedFundsSofr(kTrade).id("3M").tenor(Tenor::months(3)).build();
    CHECK(m3.legs()[0].cashflows().size() == 1);
}

TEST_CASE("DepositBuilder, FraBuilder, FutureBuilder", "[instruments][builder]") {
    Deposit d = DepositBuilder::usdSofr(kTrade).id("D-1M").tenor(Tenor::months(1)).rate(0.0383).build();
    CHECK(d.start() == Date(2026, 9, 9));
    CHECK(d.end() == Date(2026, 10, 9));
    CHECK_THAT(d.tau(), WithinAbs(30.0 / 360.0, 1e-16));
    CHECK(d.notional() == 1e6);
    Deposit d2 = DepositBuilder::usdSofr(kTrade).id("D-x").start(Date(2026, 9, 30)).tenor(Tenor::months(1)).rate(0.038).build();
    CHECK(d2.end() == Date(2026, 10, 30)); // 30 Sep is month end; EOM rule → last business day of October (Fri 30)

    Fra f = FraBuilder(rl::core::Index::termSofr3m()).id("FRA").start(Date(2026, 12, 9)).end(Date(2027, 3, 9)).rate(0.036).build();
    CHECK(f.fixingDate() == Date(2026, 12, 7)); // 2 SIFMA business days before
    CHECK_THAT(f.tau(), WithinAbs(90.0 / 360.0, 1e-16));

    Future fu = FutureBuilder::usdSofr3m().id("SR3-Z6").contract(2026, 12).price(96.33).convexityAdjustmentBp(0.3).build();
    CHECK(fu.start() == Date(2026, 12, 16));
    CHECK(fu.end() == Date(2027, 3, 17));
    CHECK_THAT(fu.tau(), WithinAbs(91.0 / 360.0, 1e-16));
    CHECK_THAT(fu.convexityAdjustment(), WithinAbs(0.3e-4, 1e-20));
    CHECK_THROWS_AS(FutureBuilder::usdSofr3m().id("x").contract(2026, 11), rl::InputError);
    Future fu2 = FutureBuilder::usdSofr3m().id("y").start(Date(2026, 9, 16)).price(96.2).build();
    CHECK(fu2.end() == Date(2026, 12, 16));
}

TEST_CASE("Leg factories: tau is under the index day count, term fixing dates respect the lag", "[instruments][builder]") {
    rl::core::ScheduleSpec spec;
    spec.effective = Date(2026, 9, 9);
    spec.termination = Date(2027, 9, 9);
    spec.frequency = Tenor::months(3);
    spec.calendar = rl::core::Calendar::get("USD-SIFMA");
    spec.paymentLag = 0;
    auto s = rl::core::Schedule::generate(spec);
    Leg t = makeTermLeg(s, rl::core::Index::termSofr3m(), 0.0, 1e6, Side::Receive, rl::core::Currency::USD);
    REQUIRE(t.cashflows().size() == 4);
    const auto& c0 = as<TermCoupon>(t.cashflows()[0]);
    CHECK(c0.fixingDate == Date(2026, 9, 4)); // 2 business days before 9 Sep, skipping Labor Day
    CHECK_THAT(c0.tau, WithinAbs(91.0 / 360.0, 1e-16));
    CHECK_THROWS_AS(makeTermLeg(s, rl::core::Index::sofr(), 0.0, 1e6, Side::Receive, rl::core::Currency::USD), rl::InputError);
    CHECK_THROWS_AS(makeOvernightLeg(s, rl::core::Index::termSofr3m(), 0.0, 1e6, Side::Receive, rl::core::Currency::USD), rl::InputError);
    CHECK_THROWS_AS(makeFixedLeg(s, rl::core::DayCount::Act360, 0.03, -1.0, Side::Pay, rl::core::Currency::USD), rl::InputError);
}
