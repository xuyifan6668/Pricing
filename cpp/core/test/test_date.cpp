#include "rl/core/date.hpp"

#include "rl/core/error.hpp"

#include <catch2/catch_test_macros.hpp>

#include <sstream>

using namespace rl::core;

// Goldens computed independently in Python: (date - date(1899,12,30)).days and strftime("%A").

TEST_CASE("Serial goldens against the Excel epoch", "[core][date]") {
    CHECK(Date(1899, 12, 30).serial() == 0);
    CHECK(Date(1899, 12, 31).serial() == 1);
    CHECK(Date(1900, 3, 1).serial() == 61); // Excel agrees from here on
    CHECK(Date(1970, 1, 1).serial() == Date::kUnixEpochSerial);
    CHECK(Date(2000, 1, 1).serial() == 36526);
    CHECK(Date(2000, 2, 29).serial() == 36585);
    CHECK(Date(2024, 1, 1).serial() == 45292);
    CHECK(Date(2026, 9, 4).serial() == 46269);
    CHECK(Date(2026, 9, 8).serial() == 46273);
    CHECK(Date(2100, 3, 1).serial() == 73110);
    CHECK(Date(1, 1, 1).serial() == Date::minDate().serial());
    CHECK(Date(9999, 12, 31).serial() == Date::maxDate().serial());
}

TEST_CASE("civilFromDays is the inverse of daysFromCivil", "[core][date]") {
    for (int s = Date::minDate().serial(); s <= Date::maxDate().serial(); s += 7919) { // prime stride
        Ymd c = Date::civilFromDays(s);
        CHECK(Date::daysFromCivil(c.year, c.month, c.day) == s);
    }
    // And densely across the region that matters.
    for (int s = Date(1990, 1, 1).serial(); s <= Date(2110, 12, 31).serial(); ++s) {
        Ymd c = Date::civilFromDays(s);
        REQUIRE(Date::daysFromCivil(c.year, c.month, c.day) == s);
    }
}

TEST_CASE("ymd round-trip and accessors", "[core][date]") {
    Date d(2026, 9, 8);
    CHECK(d.year() == 2026);
    CHECK(d.month() == 9);
    CHECK(d.day() == 8);
    CHECK(d.ymd() == Ymd{2026, 9, 8});
    CHECK(Date(2026, Month::September, 8) == d);
}

TEST_CASE("Leap years", "[core][date]") {
    CHECK(Date::isLeap(2000));
    CHECK(Date::isLeap(2024));
    CHECK_FALSE(Date::isLeap(1900));
    CHECK_FALSE(Date::isLeap(2100));
    CHECK_FALSE(Date::isLeap(2023));
    CHECK(Date::isLeap(2028));
}

TEST_CASE("Days in month sum to the year length", "[core][date]") {
    auto sum = [](int y) {
        unsigned t = 0;
        for (unsigned m = 1; m <= 12; ++m) t += Date::daysInMonth(y, m);
        return t;
    };
    CHECK(sum(2023) == 365);
    CHECK(sum(2024) == 366);
    CHECK(Date::daysInMonth(2024, 2) == 29);
    CHECK(Date::daysInMonth(2023, 2) == 28);
    CHECK(Date::daysInMonth(2026, 4) == 30);
    CHECK(Date::daysInMonth(2026, 12) == 31);
}

TEST_CASE("Weekdays", "[core][date]") {
    CHECK(Date(1899, 12, 30).weekday() == Weekday::Saturday);
    CHECK(Date(1899, 12, 31).weekday() == Weekday::Sunday);
    CHECK(Date(1970, 1, 1).weekday() == Weekday::Thursday);
    CHECK(Date(2000, 1, 1).weekday() == Weekday::Saturday);
    CHECK(Date(2026, 9, 4).weekday() == Weekday::Friday);
    CHECK(Date(2026, 9, 8).weekday() == Weekday::Tuesday);
    CHECK(Date(1, 1, 1).weekday() == Weekday::Monday); // negative serial: modulus must be non-negative
    CHECK(Date(2026, 9, 5).isWeekend());
    CHECK_FALSE(Date(2026, 9, 4).isWeekend());
}

TEST_CASE("Weekday cycles with period 7", "[core][date]") {
    int base = Date(2026, 9, 4).serial();
    for (int k = -30; k < 30; ++k) CHECK(Date(base + 7 * k).weekday() == Weekday::Friday);
}

TEST_CASE("End of month", "[core][date]") {
    CHECK(Date(2026, 2, 10).endOfMonth() == Date(2026, 2, 28));
    CHECK(Date(2028, 2, 1).endOfMonth() == Date(2028, 2, 29));
    CHECK(Date(2026, 12, 31).endOfMonth() == Date(2026, 12, 31));
    CHECK(Date(2026, 12, 31).isEndOfMonth());
    CHECK_FALSE(Date(2026, 12, 30).isEndOfMonth());
}

TEST_CASE("addMonths clamps the day and never rolls over", "[core][date]") {
    CHECK(Date(2027, 1, 31).addMonths(1) == Date(2027, 2, 28));
    CHECK(Date(2028, 1, 31).addMonths(1) == Date(2028, 2, 29));
    CHECK(Date(2027, 3, 31).addMonths(-1) == Date(2027, 2, 28));
    CHECK(Date(2026, 8, 31).addMonths(6) == Date(2027, 2, 28));
    CHECK(Date(2026, 9, 8).addMonths(12) == Date(2027, 9, 8));
    CHECK(Date(2026, 9, 8).addMonths(-9) == Date(2025, 12, 8));
    CHECK(Date(2026, 9, 8).addMonths(-21) == Date(2024, 12, 8));
    CHECK(Date(2026, 1, 15).addMonths(-1) == Date(2025, 12, 15)); // year boundary going back
    CHECK(Date(2026, 9, 8).addYears(10) == Date(2036, 9, 8));
    CHECK(Date(2028, 2, 29).addYears(1) == Date(2029, 2, 28));
}

TEST_CASE("addMonths is generated from the anchor, not iteratively", "[core][date]") {
    Date anchor(2027, 1, 31);
    CHECK(anchor.addMonths(2) == Date(2027, 3, 31)); // iterative would give 28 March
    CHECK(anchor.addMonths(3) == Date(2027, 4, 30));
}

TEST_CASE("Day arithmetic and comparisons", "[core][date]") {
    Date d(2026, 9, 8);
    CHECK(d + 1 == Date(2026, 9, 9));
    CHECK(d - 8 == Date(2026, 8, 31));
    CHECK((d + 30) - d == 30);
    CHECK(d < d + 1);
    CHECK(d >= d);
    Date e = d;
    ++e;
    e += 6;
    CHECK(e - d == 7);
}

TEST_CASE("ISO parsing and formatting", "[core][date]") {
    CHECK(Date::fromIso("2026-09-08") == Date(2026, 9, 8));
    CHECK(Date(2026, 9, 8).toIso() == "2026-09-08");
    CHECK(Date(2000, 2, 29).toIso() == "2000-02-29");
    std::ostringstream os;
    os << Date(2026, 1, 5);
    CHECK(os.str() == "2026-01-05");
    CHECK_THROWS_AS(Date::fromIso("2026/09/08"), rl::InputError);
    CHECK_THROWS_AS(Date::fromIso("20260908"), rl::InputError);
    CHECK_THROWS_AS(Date::fromIso("2026-13-01"), rl::InputError);
    CHECK_THROWS_AS(Date::fromIso("2026-02-30"), rl::InputError);
}

TEST_CASE("Constructor validation", "[core][date]") {
    CHECK_THROWS_AS(Date(2026, 0, 1), rl::InputError);
    CHECK_THROWS_AS(Date(2026, 13, 1), rl::InputError);
    CHECK_THROWS_AS(Date(2026, 2, 29), rl::InputError);
    CHECK_NOTHROW(Date(2028, 2, 29));
    CHECK_THROWS_AS(Date(2026, 4, 31), rl::InputError);
}

TEST_CASE("Enum names", "[core][date]") {
    CHECK(toString(Weekday::Monday) == "Monday");
    CHECK(toString(Month::September) == "September");
}
