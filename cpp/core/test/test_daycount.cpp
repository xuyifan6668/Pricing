#include "rl/core/daycount.hpp"

#include "rl/core/error.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using namespace rl::core;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

TEST_CASE("Day count string round trip", "[core][daycount]") {
    for (auto dc : {DayCount::Act360, DayCount::Act365Fixed, DayCount::ActActIsda, DayCount::Thirty360US, DayCount::ThirtyE360})
        CHECK(dayCountFromString(toString(dc)) == dc);
    CHECK(dayCountFromString("Act/365") == DayCount::Act365Fixed);
    CHECK(dayCountFromString("30/360 US") == DayCount::Thirty360US);
    CHECK(dayCountFromString("Eurobond") == DayCount::ThirtyE360);
    CHECK_THROWS_AS(dayCountFromString("ACT/366"), rl::InputError);
}

TEST_CASE("ACT/360 and ACT/365F are plain day ratios", "[core][daycount]") {
    Date a(2026, 9, 8), b(2027, 9, 8); // 365 actual days
    CHECK(dayCountDays(DayCount::Act360, a, b) == 365);
    CHECK_THAT(yearFraction(DayCount::Act360, a, b), WithinRel(365.0 / 360.0, 1e-15));
    CHECK_THAT(yearFraction(DayCount::Act365Fixed, a, b), WithinRel(1.0, 1e-15));
    CHECK_THAT(yearFraction(DayCount::Act365Fixed, Date(2027, 9, 8), Date(2028, 9, 8)), WithinRel(366.0 / 365.0, 1e-15));
    CHECK(yearFraction(DayCount::Act360, b, a) < 0.0);
}

TEST_CASE("ACT/ACT ISDA golden and properties", "[core][daycount]") {
    // ISDA / EMU conventions memo example: 61/365 + 121/366.
    CHECK_THAT(yearFraction(DayCount::ActActIsda, Date(2003, 11, 1), Date(2004, 5, 1)), WithinAbs(0.49772438056740775, 1e-15));
    // Same day one year later is exactly 1 unless a Feb-29 falls inside the *second* year's part.
    CHECK_THAT(yearFraction(DayCount::ActActIsda, Date(2026, 9, 8), Date(2027, 9, 8)), WithinAbs(1.0, 1e-15));
    CHECK_THAT(yearFraction(DayCount::ActActIsda, Date(2027, 9, 8), Date(2028, 9, 8)), WithinAbs(1.0008608428774608, 1e-15));
    // Within one leap year: days/366.
    CHECK_THAT(yearFraction(DayCount::ActActIsda, Date(2028, 1, 1), Date(2028, 7, 1)), WithinAbs(182.0 / 366.0, 1e-15));
    // Antisymmetry.
    CHECK_THAT(yearFraction(DayCount::ActActIsda, Date(2028, 7, 1), Date(2028, 1, 1)), WithinAbs(-182.0 / 366.0, 1e-15));
    // Telescoping: sum over consecutive periods equals the whole (exact for ACT/ACT, not ACT/365F).
    Date d0(2026, 9, 8), d1(2027, 12, 31), d2(2028, 6, 15), d3(2029, 9, 8);
    double whole = yearFraction(DayCount::ActActIsda, d0, d3);
    double parts = yearFraction(DayCount::ActActIsda, d0, d1) + yearFraction(DayCount::ActActIsda, d1, d2) + yearFraction(DayCount::ActActIsda, d2, d3);
    CHECK_THAT(parts, WithinAbs(whole, 1e-14));
}

TEST_CASE("30/360 Bond Basis rules (ISDA 2006 4.16(f))", "[core][daycount]") {
    CHECK(dayCountDays(DayCount::Thirty360US, Date(2007, 1, 31), Date(2007, 2, 28)) == 28); // D1 31→30, D2 28
    CHECK(dayCountDays(DayCount::Thirty360US, Date(2007, 1, 30), Date(2007, 3, 31)) == 60); // D2 31→30 because D1 is 30
    CHECK(dayCountDays(DayCount::Thirty360US, Date(2007, 1, 15), Date(2007, 3, 31)) == 76); // D2 stays 31: D1 < 30
    CHECK(dayCountDays(DayCount::Thirty360US, Date(2026, 9, 8), Date(2027, 9, 8)) == 360);
    CHECK(dayCountDays(DayCount::Thirty360US, Date(2026, 2, 28), Date(2026, 3, 28)) == 30); // Feb is not special
    CHECK(dayCountDays(DayCount::Thirty360US, Date(2026, 2, 28), Date(2026, 3, 31)) == 33);
    CHECK_THAT(yearFraction(DayCount::Thirty360US, Date(2026, 9, 8), Date(2027, 3, 8)), WithinRel(0.5, 1e-15));
}

TEST_CASE("30E/360 Eurobond rules (ISDA 2006 4.16(g))", "[core][daycount]") {
    CHECK(dayCountDays(DayCount::ThirtyE360, Date(2007, 1, 31), Date(2007, 2, 28)) == 28);
    CHECK(dayCountDays(DayCount::ThirtyE360, Date(2007, 1, 30), Date(2007, 3, 31)) == 60);
    CHECK(dayCountDays(DayCount::ThirtyE360, Date(2007, 1, 15), Date(2007, 3, 31)) == 75); // differs from Bond Basis
    // Any period ending on the 31st equals the same period ending on the 30th.
    for (unsigned m : {1u, 3u, 5u, 7u, 8u, 10u, 12u})
        CHECK(dayCountDays(DayCount::ThirtyE360, Date(2026, 1, 15), Date(2026, m, 31)) == dayCountDays(DayCount::ThirtyE360, Date(2026, 1, 15), Date(2026, m, 30)));
}

TEST_CASE("A full month between equal days-of-month is 30 under both 30/360 variants", "[core][daycount]") {
    for (unsigned m = 1; m <= 11; ++m) {
        CHECK(dayCountDays(DayCount::Thirty360US, Date(2026, m, 15), Date(2026, m + 1, 15)) == 30);
        CHECK(dayCountDays(DayCount::ThirtyE360, Date(2026, m, 15), Date(2026, m + 1, 15)) == 30);
    }
}
