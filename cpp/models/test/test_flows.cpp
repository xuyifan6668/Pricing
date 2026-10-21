// The bond decomposition of a swap must reprice the swap on today's curve exactly: it is the bridge
// between the trade as booked (linear pricing off the curve) and the trade as the model sees it.
#include "rl/models/flows.hpp"

#include "rl/core/error.hpp"
#include "rl/linear/annuity.hpp"
#include "rl/linear/leg_pv.hpp"
#include "rl/linear/par_rate.hpp"
#include "test_support.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using namespace rl::models;
using namespace rl::models::test;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

TEST_CASE("Swap as bonds reprices the swap off the curve to 1e-13 of notional", "[models][flows]") {
    const UsdMarket& m = usdMarket();
    for (int years : {2, 5, 10, 30}) {
        for (bool payer : {true, false}) {
            const auto swap = spotSwap(years, 0.031, payer);
            double v = 0.0;
            for (const BondFlow& f : swapAsBonds(swap, m.curves, m.asOf)) v += f.amount * f.df;
            CHECK_THAT(v, WithinAbs(rl::linear::pv(swap, m.curves, m.fixings), 1e-13));
        }
    }
}

TEST_CASE("Exercise right: forward rate and annuity equal those of the forward-starting swap", "[models][flows]") {
    const UsdMarket& m = usdMarket();
    const auto swap = spotSwap(10, 0.035);
    const auto dates = coterminalExerciseDates(swap, 1.0, 0, "USD-SIFMA");
    REQUIRE(dates.size() == 9);
    for (std::size_t i = 0; i < dates.size(); ++i) {
        const ExerciseRight r = makeExerciseRight(swap, m.curves, dates[i]);
        CHECK(r.payer);
        CHECK(r.strike == 0.035);
        CHECK(r.annuity.size() == 9 - i);
        CHECK_THAT(r.tenorYears, WithinAbs(9.0 - static_cast<double>(i), 0.02));
        // The same remaining swap built directly as a forward starter from the swap's own schedule.
        const auto& lastFixed = std::get<rl::instruments::FixedCoupon>(swap.fixedLeg().cashflows().back());
        const auto fwd = rl::instruments::OisSwapBuilder::usdSofr(m.asOf).id("fwd").effective(dates[i]).maturity(lastFixed.accrualEnd).fixedRate(0.035).notional(1.0).build();
        CHECK_THAT(r.forwardRate(), WithinAbs(rl::linear::parRate(fwd, m.curves, m.fixings), 2e-6)); // schedules may differ by a roll day
        CHECK_THAT(r.forwardValue(), WithinAbs(r.annuityValue() * (r.forwardRate() - r.strike), 1e-15));
    }
}

TEST_CASE("Payer and receiver rights are mirror images; notice days move the exercise date earlier", "[models][flows]") {
    const auto pay = rights10nc1(0.035, true), rec = rights10nc1(0.035, false);
    REQUIRE(pay.size() == rec.size());
    for (std::size_t i = 0; i < pay.size(); ++i) {
        CHECK_THAT(pay[i].forwardValue(), WithinAbs(-rec[i].forwardValue(), 1e-15));
        CHECK_THAT(pay[i].forwardRate(), WithinAbs(rec[i].forwardRate(), 1e-15));
    }
    const auto notice = rights10nc1(0.035, true, 2);
    for (std::size_t i = 0; i < pay.size(); ++i) {
        CHECK(notice[i].date < pay[i].date);
        CHECK_THAT(notice[i].forwardValue(), WithinAbs(pay[i].forwardValue(), 1e-15)); // same swap, earlier decision
    }
}

TEST_CASE("Flows validation", "[models][flows]") {
    const UsdMarket& m = usdMarket();
    const auto swap = spotSwap(5, 0.03);
    CHECK_THROWS_AS(makeExerciseRight(swap, m.curves, m.asOf), rl::InputError);                       // not in the future
    CHECK_THROWS_AS(makeExerciseRight(swap, m.curves, m.asOf.addYears(6)), rl::InputError);           // nothing left to enter
    const auto ff = rl::instruments::OisSwapBuilder::usdFedFunds(m.asOf).id("ff").tenor(rl::core::Tenor::years(5)).fixedRate(0.03).build();
    CHECK_THROWS_AS(swapAsBonds(ff, m.curves, m.asOf), rl::InputError);                               // not the discount index
}
