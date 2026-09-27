#include "rl/instruments/conventions.hpp"

#include "rl/core/error.hpp"

#include <catch2/catch_test_macros.hpp>

using namespace rl::instruments;
using rl::core::Date;
using rl::core::DayCount;
using rl::core::Tenor;

TEST_CASE("USD conventions are the market-standard ones", "[instruments][conventions]") {
    auto s = Conventions::usdSofrOis();
    CHECK(s.index.name == "USD-SOFR");
    CHECK(s.fixedDayCount == DayCount::Act360);
    CHECK(s.floatDayCount == DayCount::Act360);
    CHECK(s.fixedFrequency == Tenor::years(1));
    CHECK(s.paymentLag == 2);
    CHECK(s.spotLag == 2);
    CHECK(s.singlePeriodUpTo == Tenor::years(1));
    auto b = Conventions::usdFedFundsSofrBasis();
    CHECK(b.spreadIndex.name == "USD-FEDFUNDS");
    CHECK(b.flatIndex.name == "USD-SOFR");
    CHECK(b.spreadIndex.compounding == rl::core::Compounding::AveragedDaily);
    CHECK(b.flatIndex.compounding == rl::core::Compounding::CompoundedDaily);
    CHECK(b.frequency == Tenor::months(3));
    CHECK(Conventions::oisSwap("USD-FEDFUNDS").index.name == "USD-FEDFUNDS");
    CHECK_THROWS_AS(Conventions::oisSwap("EUR-ESTR"), rl::InputError);
    CHECK_THROWS_AS(Conventions::basisSwap("USD-SOFR", "USD-FEDFUNDS"), rl::InputError); // spread side matters
    CHECK(Conventions::future("USD-SOFR").referenceTenor == Tenor::months(3));
    CHECK_THROWS_AS(Conventions::future("USD-FEDFUNDS"), rl::InputError);
}

TEST_CASE("IMM dates: third Wednesdays, quarterly cycle, codes", "[instruments][imm]") {
    // Hand-checked against the CME calendar: Sep 2026 → 16th, Dec 2026 → 16th, Mar 2027 → 17th,
    // Jun 2027 → 16th, Sep 2027 → 15th, Jun 2028 → 21st, Sep 2028 → 20th.
    CHECK(thirdWednesday(2026, 9) == Date(2026, 9, 16));
    CHECK(thirdWednesday(2026, 12) == Date(2026, 12, 16));
    CHECK(thirdWednesday(2027, 3) == Date(2027, 3, 17));
    CHECK(thirdWednesday(2027, 6) == Date(2027, 6, 16));
    CHECK(thirdWednesday(2027, 9) == Date(2027, 9, 15));
    CHECK(thirdWednesday(2028, 6) == Date(2028, 6, 21));
    CHECK(thirdWednesday(2028, 9) == Date(2028, 9, 20));
    for (int y = 2020; y < 2040; ++y)
        for (unsigned m = 1; m <= 12; ++m) {
            Date d = thirdWednesday(y, m);
            CHECK(d.weekday() == rl::core::Weekday::Wednesday);
            CHECK(d.day() >= 15);
            CHECK(d.day() <= 21);
            CHECK(d.month() == m);
        }
    CHECK(isQuarterlyImm(Date(2026, 12, 16)));
    CHECK_FALSE(isQuarterlyImm(Date(2026, 11, 18))); // third Wednesday of November is not quarterly
    CHECK_FALSE(isQuarterlyImm(Date(2026, 12, 9)));
    CHECK(nextQuarterlyImm(Date(2026, 9, 4)) == Date(2026, 9, 16));
    CHECK(nextQuarterlyImm(Date(2026, 9, 16)) == Date(2026, 12, 16)); // strictly after
    CHECK(nextQuarterlyImm(Date(2026, 12, 20)) == Date(2027, 3, 17));
    CHECK(nthQuarterlyImm(Date(2026, 9, 4), 1) == Date(2026, 9, 16));
    CHECK(nthQuarterlyImm(Date(2026, 9, 4), 8) == Date(2028, 6, 21));
    CHECK_THROWS_AS(nthQuarterlyImm(Date(2026, 9, 4), 0), rl::InputError);
    CHECK(immMonthCode(3) == 'H');
    CHECK(immMonthCode(12) == 'Z');
    CHECK_THROWS_AS(immMonthCode(4), rl::InputError);
    CHECK(immMonthFromCode('u') == 9);
    CHECK_THROWS_AS(immMonthFromCode('X'), rl::InputError);
}

TEST_CASE("addTenor and spotDate", "[instruments][conventions]") {
    CHECK(addTenor(Date(2026, 9, 8), Tenor::weeks(1)) == Date(2026, 9, 15));
    CHECK(addTenor(Date(2026, 9, 8), Tenor::days(3)) == Date(2026, 9, 11));
    CHECK(addTenor(Date(2026, 9, 8), Tenor::months(18)) == Date(2028, 3, 8));
    CHECK(addTenor(Date(2026, 9, 8), Tenor::years(10)) == Date(2036, 9, 8));
    CHECK(addTenor(Date(2027, 1, 31), Tenor::months(1)) == Date(2027, 2, 28));
    // 2026-09-04 is a Friday; T+2 on Fedwire (Labor Day 7 Sep is a holiday) is Wednesday 9 Sep.
    CHECK(spotDate(Date(2026, 9, 4), 2, rl::core::Calendar::get("USD-FEDWIRE")) == Date(2026, 9, 9));
    CHECK(spotDate(Date(2026, 9, 4), 0, rl::core::Calendar::get("USD-FEDWIRE")) == Date(2026, 9, 4));
}
