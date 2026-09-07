#include "rl/core/date.hpp"

#include "rl/core/error.hpp"

#include <algorithm>
#include <charconv>
#include <cstdio>
#include <ostream>

namespace rl::core {

// ------------------------------------------------------------------------------------------------
// Plumbing
// ------------------------------------------------------------------------------------------------

Date::Date(int year, unsigned month, unsigned day) {
    RL_REQUIRE(month >= 1 && month <= 12, "month out of range");
    RL_REQUIRE(day >= 1 && day <= daysInMonth(year, month), "day out of range for month");
    RL_REQUIRE(year >= 1 && year <= 9999, "year out of range");
    serial_ = daysFromCivil(year, month, day);
}

Date Date::fromIso(std::string_view s) {
    RL_REQUIRE(s.size() == 10 && s[4] == '-' && s[7] == '-', "date must be YYYY-MM-DD: " + std::string(s));
    auto parse = [&](std::size_t pos, std::size_t len) {
        int v = 0;
        auto first = s.data() + pos;
        auto [p, ec] = std::from_chars(first, first + len, v);
        RL_REQUIRE(ec == std::errc{} && p == first + len, "date must be YYYY-MM-DD: " + std::string(s));
        return v;
    };
    int y = parse(0, 4);
    int m = parse(5, 2);
    int d = parse(8, 2);
    RL_REQUIRE(m >= 1 && d >= 1, "date must be YYYY-MM-DD: " + std::string(s));
    return Date(y, static_cast<unsigned>(m), static_cast<unsigned>(d));
}

std::string Date::toIso() const {
    auto [y, m, d] = ymd();
    char buf[16];
    std::snprintf(buf, sizeof buf, "%04d-%02u-%02u", y, m, d);
    return buf;
}

std::ostream& operator<<(std::ostream& os, Date d) { return os << d.toIso(); }

std::string_view toString(Weekday w) {
    static constexpr std::string_view names[] = {"Sunday",   "Monday", "Tuesday", "Wednesday",
                                                 "Thursday", "Friday", "Saturday"};
    return names[static_cast<int>(w)];
}

std::string_view toString(Month m) {
    static constexpr std::string_view names[] = {"January", "February", "March",     "April",   "May",      "June",
                                                 "July",    "August",   "September", "October", "November", "December"};
    return names[static_cast<unsigned>(m) - 1];
}

// ------------------------------------------------------------------------------------------------
// Calendar arithmetic (Hinnant date algorithms)
// ------------------------------------------------------------------------------------------------

/**
 * Gregorian leap-year rule: divisible by 4, except centuries, except those divisible by 400.
 *   leap(y) = (y mod 4 = 0) ∧ (y mod 100 ≠ 0 ∨ y mod 400 = 0)
 * 1900 is not leap, 2000 is, 2100 is not. Called from the schedule builder and day counts, so keep it
 * branch-light. Test: the three years above plus 2024 (leap), 2023 (not).
 */
bool Date::isLeap(int year) { return (year % 4 == 0) && ((year % 100 != 0) || (year % 400 == 0)); }

/**
 * Days in month m of year y: 31 for {1,3,5,7,8,10,12}, 30 for {4,6,9,11}, 28/29 for February.
 * Precondition: 1 ≤ m ≤ 12 (unchecked here; the Date ctor checks). Test: every month of 2023 and
 * 2024 sums to 365 / 366.
 */
unsigned Date::daysInMonth(int year, unsigned month) {
    static constexpr unsigned daysPerMonth[] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
    return daysPerMonth[month - 1] + ((month == 2 && isLeap(year)) ? 1u : 0u);
}

namespace {
// C++ `/` and `%` truncate toward zero; the date algorithms need floor semantics for negative inputs.
constexpr int floorDiv(int a, int b) {
    const int q = a / b, r = a % b;
    return (r != 0 && ((r < 0) != (b < 0))) ? q - 1 : q;
}
constexpr int floorMod(int a, int b) {
    const int r = a % b;
    return r < 0 ? r + b : r;
}
} // namespace
/**
 * Proleptic-Gregorian civil date -> serial day, epoch 1899-12-30.
 *
 * Use Howard Hinnant's `days_from_civil` (http://howardhinnant.github.io/date_algorithms.html,
 * "days_from_civil"), which gives days since 1970-01-01 with no loops and no tables:
 *
 *   y' = y − [m ≤ 2]                      (shift the year so March is month 0; Feb 29 is then last)
 *   era = ⌊y'/400⌋                        (floor division — careful with negative y')
 *   yoe = y' − 400·era                     ∈ [0, 399]
 *   doy = ⌊(153·(m + (m > 2 ? −3 : 9)) + 2)/5⌋ + d − 1     ∈ [0, 365]
 *   doe = yoe·365 + ⌊yoe/4⌋ − ⌊yoe/100⌋ + doy               ∈ [0, 146096]
 *   days_since_1970 = era·146097 + doe − 719468
 *
 * Then add kUnixEpochSerial (25569) to land on our epoch. The 153/5 trick encodes the 31,30,31,30,31
 * month-length pattern starting from March; 146097 is the number of days in a 400-year cycle.
 *
 * Numerical: all integer; watch the floor of negative numbers (C++ `/` truncates toward zero, so
 * handle y' < 0 explicitly). Test: round-trip with civilFromDays for every serial in [-5000, 3e6],
 * and the goldens in test_date.cpp (1899-12-30 → 0, 1970-01-01 → 25569, 2000-01-01 → 36526).
 */
int Date::daysFromCivil(int year, unsigned month, unsigned day) {
    const int y = year - (month <= 2 ? 1 : 0);
    const int era = floorDiv(y, 400);
    const int yoe = y - 400 * era;
    const int mp = static_cast<int>(month) + (month > 2 ? -3 : 9);
    const int doy = (153 * mp + 2) / 5 + static_cast<int>(day) - 1;
    const int doe = yoe * 365 + yoe / 4 - yoe / 100 + doy;
    return era * 146097 + doe - 719468 + kUnixEpochSerial;
}

/**
 * Serial day -> civil date. Hinnant's `civil_from_days`, the exact inverse of the above:
 *
 *   z   = serial − 25569 + 719468            (days since 0000-03-01)
 *   era = ⌊z/146097⌋
 *   doe = z − era·146097                     ∈ [0, 146096]
 *   yoe = ⌊(doe − ⌊doe/1460⌋ + ⌊doe/36524⌋ − ⌊doe/146096⌋)/365⌋   ∈ [0, 399]
 *   y   = yoe + era·400
 *   doy = doe − (365·yoe + ⌊yoe/4⌋ − ⌊yoe/100⌋)                  ∈ [0, 365]
 *   mp  = ⌊(5·doy + 2)/153⌋                                       ∈ [0, 11], 0 = March
 *   d   = doy − ⌊(153·mp + 2)/5⌋ + 1
 *   m   = mp < 10 ? mp + 3 : mp − 9
 *   return {y + (m ≤ 2), m, d}
 *
 * Hot path: called by every day-count and schedule step; O(1), no allocation.
 */
Ymd Date::civilFromDays(int serial) {
    const int z = serial - kUnixEpochSerial + 719468;
    const int era = floorDiv(z, 146097);
    const int doe = z - era * 146097;                                      // [0, 146096]
    const int yoe = (doe - doe / 1460 + doe / 36524 - doe / 146096) / 365; // [0, 399]
    const int y = yoe + era * 400;
    const int doy = doe - (365 * yoe + yoe / 4 - yoe / 100);               // [0, 365]
    const int mp = (5 * doy + 2) / 153;                                    // [0, 11], 0 = March
    const auto d = static_cast<unsigned>(doy - (153 * mp + 2) / 5 + 1);
    const auto m = static_cast<unsigned>(mp < 10 ? mp + 3 : mp - 9);
    return Ymd{y + (m <= 2 ? 1 : 0), m, d};
}
/**
 * Day of week. Serial 0 (1899-12-30) is a Saturday, so with Sunday = 0:
 *   weekday = (serial + 6) mod 7,  using a *non-negative* modulus for serial < 0.
 * Test: 2026-09-04 is a Friday; 1970-01-01 is a Thursday; 1899-12-31 is a Sunday.
 */
Weekday Date::weekday() const { return static_cast<Weekday>(floorMod(serial_ + 6, 7)); }

/**
 * Last calendar day of this date's month: {y, m, daysInMonth(y, m)}.
 */
Date Date::endOfMonth() const {
    const Ymd c = civilFromDays(serial_);
    return Date(c.year, c.month, daysInMonth(c.year, c.month));
}

/**
 * Add n months (n may be negative), clamping the day to the target month's length:
 *   2027-01-31 + 1M = 2027-02-28,  2028-01-31 + 1M = 2028-02-29,  2027-03-31 − 1M = 2027-02-28.
 * Year/month arithmetic: total = 12·y + (m − 1) + n; y' = ⌊total/12⌋; m' = total mod 12 + 1 — again
 * with floor semantics for negatives. The day is clamped, NOT rolled into the next month.
 *
 * Schedule generation must call this from the *anchor* date each time (anchor + k·n), never
 * iteratively from the previous clamped date; otherwise 31-Jan → 28-Feb → 28-Mar loses the 31.
 * Test: the three cases above and 2026-08-31 + 6M = 2027-02-28.
 */
Date Date::addMonths(int n) const {
    const Ymd c = civilFromDays(serial_);
    const int total = static_cast<int>(c.month) - 1 + n;
    const int year = c.year + floorDiv(total, 12);
    const auto month = static_cast<unsigned>(floorMod(total, 12)) + 1;
    return Date(year, month, std::min(c.day, daysInMonth(year, month)));
}

} // namespace rl::core
