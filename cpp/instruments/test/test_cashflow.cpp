// The cashflow variant, Leg and the trade types.
#include "rl/instruments/basis_swap.hpp"
#include "rl/instruments/bermudan.hpp"
#include "rl/instruments/cashflow.hpp"
#include "rl/instruments/deposit.hpp"
#include "rl/instruments/fra.hpp"
#include "rl/instruments/future.hpp"
#include "rl/instruments/leg.hpp"
#include "rl/instruments/swap.hpp"
#include "rl/instruments/swaption.hpp"

#include "rl/core/error.hpp"

#include <catch2/catch_test_macros.hpp>

#include <type_traits>
#include <variant>

using namespace rl::instruments;
using rl::core::Currency;
using rl::core::Date;
using rl::core::Index;

namespace {
FixedCoupon fixed(Date s, Date e, Date p, double tau, double n, double r) {
    FixedCoupon c;
    c.accrualStart = s; c.accrualEnd = e; c.payment = p; c.tau = tau; c.notional = n; c.rate = r;
    return c;
}
OvernightCoupon ois(Date s, Date e, Date p, double tau, double n, double spread = 0.0) {
    OvernightCoupon c;
    c.accrualStart = s; c.accrualEnd = e; c.payment = p; c.tau = tau; c.notional = n; c.index = Index::sofr();
    c.spread = spread; c.lookbackDays = 0; c.lockoutDays = 0;
    return c;
}
Leg twoFixed(Side side) {
    return Leg({fixed(Date(2026, 9, 9), Date(2027, 9, 9), Date(2027, 9, 13), 1.0139, 1e6, 0.035),
                fixed(Date(2027, 9, 9), Date(2028, 9, 11), Date(2028, 9, 13), 1.0194, 1e6, 0.035)}, side, Currency::USD);
}
Leg twoOis(Side side) {
    return Leg({ois(Date(2026, 9, 9), Date(2027, 9, 9), Date(2027, 9, 13), 1.0139, 1e6),
                ois(Date(2027, 9, 9), Date(2028, 9, 11), Date(2028, 9, 13), 1.0194, 1e6)}, side, Currency::USD);
}
} // namespace

TEST_CASE("Cashflow variant: closed set, visit with a sum functor", "[instruments][cashflow]") {
    static_assert(std::variant_size_v<Cashflow> == 4, "the cashflow set is closed at four alternatives in v1");
    static_assert(std::is_default_constructible_v<FixedCoupon> && std::is_default_constructible_v<OvernightCoupon> &&
                  std::is_default_constructible_v<TermCoupon> && std::is_default_constructible_v<NotionalExchange>);
    std::vector<Cashflow> cfs{fixed(Date(2026, 9, 9), Date(2027, 9, 9), Date(2027, 9, 13), 1.0, 100.0, 0.03),
                              ois(Date(2026, 9, 9), Date(2027, 9, 9), Date(2027, 9, 13), 1.0, 200.0),
                              NotionalExchange{Date(2027, 9, 13), -50.0}};
    struct Notional {
        double operator()(const FixedCoupon& c) const { return c.notional; }
        double operator()(const OvernightCoupon& c) const { return c.notional; }
        double operator()(const TermCoupon& c) const { return c.notional; }
        double operator()(const NotionalExchange& c) const { return c.amount; }
    };
    double total = 0.0;
    for (const auto& cf : cfs) total += std::visit(Notional{}, cf);
    CHECK(total == 250.0);
}

TEST_CASE("Leg: sign, side, currency, ordering validation", "[instruments][leg]") {
    Leg r = twoFixed(Side::Receive);
    CHECK(r.sign() == 1.0);
    CHECK(r.side() == Side::Receive);
    CHECK(r.currency() == Currency::USD);
    CHECK(r.cashflows().size() == 2);
    CHECK_FALSE(r.empty());
    CHECK(r.lastPayment() == Date(2028, 9, 13));
    CHECK(twoFixed(Side::Pay).sign() == -1.0);
    // Payment dates must be non-decreasing.
    CHECK_THROWS_AS(Leg({fixed(Date(2027, 9, 9), Date(2028, 9, 11), Date(2028, 9, 13), 1.0, 1e6, 0.03),
                         fixed(Date(2026, 9, 9), Date(2027, 9, 9), Date(2027, 9, 13), 1.0, 1e6, 0.03)},
                        Side::Pay, Currency::USD),
                    rl::InputError);
}

TEST_CASE("Swap: legs must be opposite-sided and one-currency; maturity is the last payment", "[instruments][swap]") {
    Swap s("S1", twoFixed(Side::Pay), twoOis(Side::Receive));
    CHECK(s.id() == "S1");
    CHECK(s.fixedLeg().sign() == -1.0);
    CHECK(s.floatLeg().sign() == 1.0);
    CHECK(s.maturity() == Date(2028, 9, 13));
    CHECK_THROWS_AS(Swap("bad", twoFixed(Side::Pay), twoOis(Side::Pay)), rl::InputError);
    Leg eur = Leg({ois(Date(2026, 9, 9), Date(2027, 9, 9), Date(2027, 9, 13), 1.0, 1e6)}, Side::Receive, Currency::EUR);
    CHECK_THROWS_AS(Swap("bad", twoFixed(Side::Pay), eur), rl::InputError);
}

TEST_CASE("BasisSwap, Deposit, Fra, Future carry what the residuals need", "[instruments][trades]") {
    BasisSwap b("B1", twoOis(Side::Pay), twoOis(Side::Receive));
    CHECK(b.id() == "B1");
    CHECK(b.legs()[0].sign() == -1.0);
    CHECK(b.legs()[1].sign() == 1.0);
    CHECK(b.maturity() == Date(2028, 9, 13));
    CHECK_THROWS_AS(BasisSwap("bad", twoOis(Side::Pay), twoOis(Side::Pay)), rl::InputError);

    Deposit d("D1", Index::sofr(), Date(2026, 9, 9), Date(2026, 10, 9), 30.0 / 360.0, 1e6, 0.0383);
    CHECK(d.id() == "D1");
    CHECK(d.index().name == "USD-SOFR");
    CHECK(d.start() == Date(2026, 9, 9));
    CHECK(d.end() == Date(2026, 10, 9));
    CHECK(d.tau() == 30.0 / 360.0);
    CHECK(d.notional() == 1e6);
    CHECK(d.rate() == 0.0383);

    Fra f("F1", Index::termSofr3m(), Date(2026, 12, 7), Date(2026, 12, 9), Date(2027, 3, 9), 90.0 / 360.0, 1e6, 0.036);
    CHECK(f.fixingDate() == Date(2026, 12, 7));
    CHECK(f.start() == Date(2026, 12, 9));
    CHECK(f.end() == Date(2027, 3, 9));
    CHECK(f.rate() == 0.036);

    Future fu("SR3-Z6", Index::sofr(), Date(2026, 12, 16), Date(2027, 3, 17), 91.0 / 360.0, 96.33, 0.3e-4);
    CHECK(fu.id() == "SR3-Z6");
    CHECK(fu.start() == Date(2026, 12, 16));
    CHECK(fu.end() == Date(2027, 3, 17));
    CHECK(fu.price() == 96.33);
    CHECK(fu.convexityAdjustment() == 0.3e-4);
    CHECK(fu.tau() == 91.0 / 360.0);
}

TEST_CASE("Swaption and BermudanSwaption own their underlying", "[instruments][swaption]") {
    Swap s("S1", twoFixed(Side::Pay), twoOis(Side::Receive));
    Swaption sw("SW1", s, Date(2026, 9, 7), Settlement::Physical);
    CHECK(sw.id() == "SW1");
    CHECK(sw.expiry() == Date(2026, 9, 7));
    CHECK(sw.settlement() == Settlement::Physical);
    CHECK(sw.underlying().id() == "S1");
    BermudanSwaption bm("BM1", s, {Date(2026, 9, 7), Date(2027, 9, 7)}, 2);
    CHECK(bm.exerciseDates().size() == 2);
    CHECK(bm.noticeDays() == 2);
    CHECK(bm.underlying().fixedLeg().cashflows().size() == 2);
    static_assert(!std::is_same_v<decltype(std::declval<const Swaption&>().underlying()), Swap>,
                  "underlying() should return a const reference, not a copy");
}
