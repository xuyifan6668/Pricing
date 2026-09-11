#include "rl/core/schedule.hpp"

#include "rl/core/error.hpp"

#include <catch2/catch_test_macros.hpp>

#include <vector>

using namespace rl::core;
using BDC = BusinessDayConvention;

namespace {
std::vector<Date> iso(std::initializer_list<const char*> xs) {
    std::vector<Date> out;
    for (auto s : xs) out.push_back(Date::fromIso(s));
    return out;
}
ScheduleSpec sofrSpec() {
    ScheduleSpec s;
    s.effective = Date(2026, 9, 8);
    s.termination = Date(2036, 9, 8);
    s.frequency = Tenor::years(1);
    s.calendar = Calendar::get("USD-SIFMA");
    s.convention = BDC::ModifiedFollowing;
    s.terminationConvention = BDC::ModifiedFollowing;
    s.stub = StubRule::ShortFront;
    s.endOfMonth = true; // no effect: the 8th is not month end
    s.paymentLag = 2;
    return s;
}
} // namespace

// The end-to-end check: a 10y USD SOFR swap from 2026-09-08, annual, MF, EOM, 2-day
// payment lag, matches a hand-built schedule date for date. Weekdays checked by hand:
// 8 Sep is Sat in 2029 and 2035 (→ Mon 10), Sun in 2030 (→ Mon 9); Labor Day never collides.
TEST_CASE("10y annual SOFR schedule from 2026-09-08", "[core][schedule]") {
    Schedule s = Schedule::generate(sofrSpec());
    REQUIRE(s.size() == 10);
    CHECK(s.unadjustedDates() == iso({"2026-09-08", "2027-09-08", "2028-09-08", "2029-09-08", "2030-09-08", "2031-09-08",
                                      "2032-09-08", "2033-09-08", "2034-09-08", "2035-09-08", "2036-09-08"}));
    CHECK(s.dates() == iso({"2026-09-08", "2027-09-08", "2028-09-08", "2029-09-10", "2030-09-09", "2031-09-08",
                            "2032-09-08", "2033-09-08", "2034-09-08", "2035-09-10", "2036-09-08"}));
    CHECK(s.paymentDates() == iso({"2027-09-10", "2028-09-12", "2029-09-12", "2030-09-11", "2031-09-10",
                                   "2032-09-10", "2033-09-12", "2034-09-12", "2035-09-12", "2036-09-10"}));
    for (const auto& p : s) CHECK(p.regular);
    CHECK(s.front().start == Date(2026, 9, 8));
    CHECK(s.back().end == Date(2036, 9, 8));
    CHECK(s[3].start == Date(2029, 9, 10));
    CHECK(s[3].end == Date(2030, 9, 9));
}

TEST_CASE("Monthly from 31 January: with EOM every date is month end, without it the 31st is kept where possible", "[core][schedule]") {
    ScheduleSpec s;
    s.effective = Date(2027, 1, 31);
    s.termination = Date(2027, 7, 31);
    s.frequency = Tenor::months(1);
    s.calendar = Calendar::weekendsOnly();
    s.convention = BDC::Unadjusted;
    s.terminationConvention = BDC::Unadjusted;
    s.stub = StubRule::ShortBack; // forward generation from the effective date
    s.endOfMonth = false;
    CHECK(Schedule::generate(s).unadjustedDates() == iso({"2027-01-31", "2027-02-28", "2027-03-31", "2027-04-30", "2027-05-31", "2027-06-30", "2027-07-31"}));
    s.endOfMonth = true;
    CHECK(Schedule::generate(s).unadjustedDates() == iso({"2027-01-31", "2027-02-28", "2027-03-31", "2027-04-30", "2027-05-31", "2027-06-30", "2027-07-31"}));
    // Backward generation (front stub rule, no actual stub) from a 31-Jul anchor gives the same dates.
    s.stub = StubRule::ShortFront;
    CHECK(Schedule::generate(s).unadjustedDates() == iso({"2027-01-31", "2027-02-28", "2027-03-31", "2027-04-30", "2027-05-31", "2027-06-30", "2027-07-31"}));
    // Anchor on 30 April going backwards monthly: 30-Mar, 28-Feb, 30-Jan (not 31) without EOM.
    s.effective = Date(2027, 1, 30);
    s.termination = Date(2027, 4, 30);
    s.endOfMonth = false;
    CHECK(Schedule::generate(s).unadjustedDates() == iso({"2027-01-30", "2027-02-28", "2027-03-30", "2027-04-30"}));
    // With EOM, 30 April is month end, so every generated date is a month end: 31-Mar, 28-Feb, 31-Jan.
    s.effective = Date(2027, 1, 31);
    s.endOfMonth = true;
    CHECK(Schedule::generate(s).unadjustedDates() == iso({"2027-01-31", "2027-02-28", "2027-03-31", "2027-04-30"}));
    // EOM with a 30-Jan effective date leaves a one-day front stub (30 → 31 Jan); that is what the
    // rule says and what every production scheduler does. Passing effective = 31-Jan is the fix.
    s.effective = Date(2027, 1, 30);
    CHECK(Schedule::generate(s).unadjustedDates() == iso({"2027-01-30", "2027-01-31", "2027-02-28", "2027-03-31", "2027-04-30"}));
    CHECK_FALSE(Schedule::generate(s).front().regular);
}

TEST_CASE("Leap-year February starts", "[core][schedule]") {
    ScheduleSpec s;
    s.effective = Date(2028, 2, 29);
    s.termination = Date(2029, 2, 28);
    s.frequency = Tenor::months(6);
    s.calendar = Calendar::weekendsOnly();
    s.convention = BDC::Unadjusted;
    s.terminationConvention = BDC::Unadjusted;
    s.stub = StubRule::ShortBack;
    s.endOfMonth = true;
    CHECK(Schedule::generate(s).unadjustedDates() == iso({"2028-02-29", "2028-08-31", "2029-02-28"}));
    s.endOfMonth = false;
    CHECK(Schedule::generate(s).unadjustedDates() == iso({"2028-02-29", "2028-08-29", "2029-02-28"}));
    s.effective = Date(2027, 2, 28);
    s.termination = Date(2028, 2, 29);
    s.stub = StubRule::ShortFront; // backward from 29-Feb-2028
    s.endOfMonth = true;
    CHECK(Schedule::generate(s).unadjustedDates() == iso({"2027-02-28", "2027-08-31", "2028-02-29"}));
}

TEST_CASE("Quarterly EOM from 30 November 2026 with MF adjustment over a Memorial Day month end", "[core][schedule]") {
    ScheduleSpec s;
    s.effective = Date(2026, 11, 30);
    s.termination = Date(2027, 11, 30);
    s.frequency = Tenor::months(3);
    s.calendar = Calendar::get("USD-SIFMA");
    s.stub = StubRule::ShortFront;
    s.endOfMonth = true;
    Schedule sch = Schedule::generate(s);
    CHECK(sch.unadjustedDates() == iso({"2026-11-30", "2027-02-28", "2027-05-31", "2027-08-31", "2027-11-30"}));
    // 28 Feb 2027 is a Sunday → MF → Fri 26 Feb. 31 May 2027 is Memorial Day → MF → Fri 28 May.
    CHECK(sch.dates() == iso({"2026-11-30", "2027-02-26", "2027-05-28", "2027-08-31", "2027-11-30"}));
    for (const auto& p : sch) CHECK(p.regular);
}

TEST_CASE("Stub rules on a 3y8m semiannual swap", "[core][schedule]") {
    ScheduleSpec s;
    s.effective = Date(2026, 9, 8);
    s.termination = Date(2030, 5, 8);
    s.frequency = Tenor::months(6);
    s.calendar = Calendar::weekendsOnly();
    s.convention = BDC::Unadjusted;
    s.terminationConvention = BDC::Unadjusted;

    SECTION("ShortFront") {
        s.stub = StubRule::ShortFront;
        Schedule sch = Schedule::generate(s);
        CHECK(sch.unadjustedDates() == iso({"2026-09-08", "2026-11-08", "2027-05-08", "2027-11-08", "2028-05-08", "2028-11-08", "2029-05-08", "2029-11-08", "2030-05-08"}));
        CHECK_FALSE(sch.front().regular);
        CHECK(sch[1].regular);
        CHECK(sch.back().regular);
    }
    SECTION("LongFront") {
        s.stub = StubRule::LongFront;
        Schedule sch = Schedule::generate(s);
        CHECK(sch.unadjustedDates() == iso({"2026-09-08", "2027-05-08", "2027-11-08", "2028-05-08", "2028-11-08", "2029-05-08", "2029-11-08", "2030-05-08"}));
        CHECK_FALSE(sch.front().regular);
    }
    SECTION("ShortBack") {
        s.stub = StubRule::ShortBack;
        Schedule sch = Schedule::generate(s);
        CHECK(sch.unadjustedDates() == iso({"2026-09-08", "2027-03-08", "2027-09-08", "2028-03-08", "2028-09-08", "2029-03-08", "2029-09-08", "2030-03-08", "2030-05-08"}));
        CHECK(sch.front().regular);
        CHECK_FALSE(sch.back().regular);
    }
    SECTION("LongBack") {
        s.stub = StubRule::LongBack;
        Schedule sch = Schedule::generate(s);
        CHECK(sch.unadjustedDates() == iso({"2026-09-08", "2027-03-08", "2027-09-08", "2028-03-08", "2028-09-08", "2029-03-08", "2029-09-08", "2030-05-08"}));
        CHECK_FALSE(sch.back().regular);
    }
    SECTION("None throws when the dates do not fit") {
        s.stub = StubRule::None;
        CHECK_THROWS_AS(Schedule::generate(s), rl::InputError);
    }
    SECTION("None works when they do") {
        s.stub = StubRule::None;
        s.termination = Date(2030, 3, 8);
        CHECK(Schedule::generate(s).size() == 7);
    }
    SECTION("No stub needed: front and back rules agree") {
        s.termination = Date(2030, 3, 8);
        s.stub = StubRule::ShortFront;
        auto a = Schedule::generate(s).unadjustedDates();
        s.stub = StubRule::LongBack;
        auto b = Schedule::generate(s).unadjustedDates();
        CHECK(a == b);
        CHECK(a.size() == 8);
        for (const auto& p : Schedule::generate(s)) CHECK(p.regular);
    }
}

TEST_CASE("Explicit first regular date and next-to-last date", "[core][schedule]") {
    ScheduleSpec s;
    s.effective = Date(2026, 9, 8);
    s.termination = Date(2030, 5, 8);
    s.frequency = Tenor::months(6);
    s.calendar = Calendar::weekendsOnly();
    s.convention = BDC::Unadjusted;
    s.terminationConvention = BDC::Unadjusted;
    s.stub = StubRule::ShortBack;
    s.firstRegularDate = Date(2026, 11, 8);
    // Front stub to 8 Nov, regular forward from there, remainder at the back.
    CHECK(Schedule::generate(s).unadjustedDates() == iso({"2026-09-08", "2026-11-08", "2027-05-08", "2027-11-08", "2028-05-08", "2028-11-08", "2029-05-08", "2029-11-08", "2030-05-08"}));
    s.firstRegularDate.reset();
    s.stub = StubRule::ShortFront;
    s.nextToLastDate = Date(2030, 3, 8);
    CHECK(Schedule::generate(s).unadjustedDates() == iso({"2026-09-08", "2027-03-08", "2027-09-08", "2028-03-08", "2028-09-08", "2029-03-08", "2029-09-08", "2030-03-08", "2030-05-08"}));
    s.nextToLastDate = Date(2031, 1, 1);
    CHECK_THROWS_AS(Schedule::generate(s), rl::InputError);
}

TEST_CASE("Roll day", "[core][schedule]") {
    ScheduleSpec s;
    s.effective = Date(2026, 9, 8);
    s.termination = Date(2027, 9, 8);
    s.frequency = Tenor::months(3);
    s.calendar = Calendar::weekendsOnly();
    s.convention = BDC::Unadjusted;
    s.terminationConvention = BDC::Unadjusted;
    s.stub = StubRule::ShortFront;
    s.rollDay = 20;
    // Generated dates sit on the 20th; effective and termination are untouched; front stub 8 Sep → 20 Sep... no:
    // backward from termination: anchors 20-Jun-2027, 20-Mar-2027, 20-Dec-2026, 20-Sep-2026, then the effective date.
    CHECK(Schedule::generate(s).unadjustedDates() == iso({"2026-09-08", "2026-09-20", "2026-12-20", "2027-03-20", "2027-06-20", "2027-09-08"}));
    s.rollDay = 31; // "last day"
    CHECK(Schedule::generate(s).unadjustedDates() == iso({"2026-09-08", "2026-09-30", "2026-12-31", "2027-03-31", "2027-06-30", "2027-09-08"}));
    s.rollDay = 0;
    CHECK_THROWS_AS(Schedule::generate(s), rl::InputError);
}

TEST_CASE("Weekly and daily frequencies add calendar days", "[core][schedule]") {
    ScheduleSpec s;
    s.effective = Date(2026, 9, 1);
    s.termination = Date(2026, 9, 29);
    s.frequency = Tenor::weeks(1);
    s.calendar = Calendar::weekendsOnly();
    s.convention = BDC::Unadjusted;
    s.terminationConvention = BDC::Unadjusted;
    s.stub = StubRule::ShortFront;
    CHECK(Schedule::generate(s).unadjustedDates() == iso({"2026-09-01", "2026-09-08", "2026-09-15", "2026-09-22", "2026-09-29"}));
}

TEST_CASE("Specification validation", "[core][schedule]") {
    ScheduleSpec s;
    s.effective = Date(2026, 9, 8);
    s.termination = Date(2026, 9, 8);
    CHECK_THROWS_AS(Schedule::generate(s), rl::InputError);
    s.termination = Date(2025, 9, 8);
    CHECK_THROWS_AS(Schedule::generate(s), rl::InputError);
    s.termination = Date(2027, 9, 8);
    s.frequency = Tenor::months(0);
    CHECK_THROWS_AS(Schedule::generate(s), rl::InputError);
}

TEST_CASE("build() plumbing: adjustment, payment lag, regular flags", "[core][schedule]") {
    ScheduleSpec s = sofrSpec();
    s.termination = Date(2028, 9, 8);
    Schedule sch = Schedule::build(s, iso({"2026-09-08", "2027-09-08", "2028-09-08"}));
    REQUIRE(sch.size() == 2);
    CHECK(sch[0].payment == Date(2027, 9, 10));
    CHECK(sch[1].payment == Date(2028, 9, 12));
    CHECK(sch[0].regular);
    CHECK_THROWS_AS(Schedule::build(s, iso({"2026-09-08", "2028-09-08", "2027-09-08"})), rl::InputError);
    CHECK_THROWS_AS(Schedule::build(s, iso({"2026-09-09", "2028-09-08"})), rl::InputError);
}

TEST_CASE("Stub rule strings", "[core][schedule]") {
    CHECK(stubRuleFromString("ShortFront") == StubRule::ShortFront);
    CHECK(stubRuleFromString("LB") == StubRule::LongBack);
    CHECK(toString(StubRule::ShortBack) == "ShortBack");
    CHECK_THROWS_AS(stubRuleFromString("Middle"), rl::InputError);
}
