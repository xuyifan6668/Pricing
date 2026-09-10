#include "rl/core/daycount.hpp"

#include "rl/core/error.hpp"

#include <algorithm>
#include <cctype>
#include <string>

namespace rl::core {

// ------------------------------------------------------------------------------------------------
// Plumbing
// ------------------------------------------------------------------------------------------------

std::string_view toString(DayCount dc) {
    switch (dc) {
    case DayCount::Act360: return "ACT/360";
    case DayCount::Act365Fixed: return "ACT/365F";
    case DayCount::ActActIsda: return "ACT/ACT";
    case DayCount::Thirty360US: return "30/360";
    case DayCount::ThirtyE360: return "30E/360";
    }
    return "?";
}

DayCount dayCountFromString(std::string_view in) {
    std::string s;
    for (char c : in)
        if (!std::isspace(static_cast<unsigned char>(c))) s.push_back(static_cast<char>(std::toupper(static_cast<unsigned char>(c))));
    if (s == "ACT/360" || s == "A/360" || s == "ACTUAL/360") return DayCount::Act360;
    if (s == "ACT/365F" || s == "ACT/365" || s == "A/365F" || s == "ACT/365FIXED" || s == "ACTUAL/365(FIXED)")
        return DayCount::Act365Fixed;
    if (s == "ACT/ACT" || s == "ACT/ACTISDA" || s == "A/A" || s == "ACTUAL/ACTUAL") return DayCount::ActActIsda;
    if (s == "30/360" || s == "30/360US" || s == "BOND" || s == "BONDBASIS" || s == "30U/360") return DayCount::Thirty360US;
    if (s == "30E/360" || s == "EUROBOND" || s == "30/360E") return DayCount::ThirtyE360;
    throw InputError("unknown day count: " + std::string(in));
}

int dayCountDays(DayCount dc, Date start, Date end) {
    switch (dc) {
    case DayCount::Act360:
    case DayCount::Act365Fixed:
    case DayCount::ActActIsda: return end - start;
    case DayCount::Thirty360US: return detail::days30360US(start, end);
    case DayCount::ThirtyE360: return detail::days30E360(start, end);
    }
    throw InputError("bad DayCount enum");
}

double yearFraction(DayCount dc, Date start, Date end) {
    switch (dc) {
    case DayCount::Act360: return static_cast<double>(end - start) / 360.0;
    case DayCount::Act365Fixed: return static_cast<double>(end - start) / 365.0;
    case DayCount::ActActIsda: return detail::yfActActIsda(start, end);
    case DayCount::Thirty360US: return static_cast<double>(detail::days30360US(start, end)) / 360.0;
    case DayCount::ThirtyE360: return static_cast<double>(detail::days30E360(start, end)) / 360.0;
    }
    throw InputError("bad DayCount enum");
}

// ------------------------------------------------------------------------------------------------
// Conventions (ISDA 2006 §4.16 transcribed, with the goldens in test_daycount.cpp)
// ------------------------------------------------------------------------------------------------

namespace detail {

/**
 * 30/360 (Bond Basis), ISDA 2006 §4.16(f).
 *
 *   days = 360·(Y₂ − Y₁) + 30·(M₂ − M₁) + (D₂' − D₁')
 *   D₁' = 30 if D₁ = 31, else D₁
 *   D₂' = 30 if D₂ = 31 and D₁' = 30 (i.e. D₁ ∈ {30, 31}), else D₂
 *
 * Note the asymmetry: D₂ = 31 is only pulled back when D₁ was already 30/31. So
 * 2007-01-15 → 2007-03-31 is 76 days here but 75 under 30E/360. February is NOT special in this
 * convention (the SIA "30/360 US" variant with end-of-February rules is a different beast:
 * // SIMPLIFICATION: not provided; USD swaps use Bond Basis, which is what this is).
 *
 * Test: goldens in test_daycount.cpp; a full month between the same day-of-month is always 30.
 */
int days30360US(Date start, Date end) {
    int y1 = start.year();
    unsigned m1 = start.month();
    unsigned d1 = start.day();

    int y2 = end.year();
    unsigned m2 = end.month();
    unsigned d2 = end.day();

    // D1' = 30 if D1 = 31, else D1
    unsigned d1p = (d1 == 31) ? 30 : d1;

    // D2' = 30 if D2 = 31 AND D1' ∈ {30, 31} (i.e., D1' == 30), else D2
    unsigned d2p = (d2 == 31 && d1p == 30) ? 30 : d2;

    return 360 * (y2 - y1) + 30 * (static_cast<int>(m2) - static_cast<int>(m1))
           + (static_cast<int>(d2p) - static_cast<int>(d1p));
}

/**
 * 30E/360 (Eurobond Basis), ISDA 2006 §4.16(g).
 *
 *   D₁' = min(D₁, 30),  D₂' = min(D₂, 30)          (both 31s become 30, unconditionally)
 *   days = 360·(Y₂ − Y₁) + 30·(M₂ − M₁) + (D₂' − D₁')
 *
 * Test: every period ending on the 31st equals the same period ending on the 30th.
 */
int days30E360(Date start, Date end) {
    int y1 = start.year();
    unsigned m1 = start.month();
    unsigned d1 = start.day();

    int y2 = end.year();
    unsigned m2 = end.month();
    unsigned d2 = end.day();

    // Both ends: 31 → 30, unconditionally (this is what distinguishes Eurobond from Bond Basis).
    unsigned d1p = (d1 > 30) ? 30 : d1;
    unsigned d2p = (d2 > 30) ? 30 : d2;

    return 360 * (y2 - y1) + 30 * (static_cast<int>(m2) - static_cast<int>(m1))
           + (static_cast<int>(d2p) - static_cast<int>(d1p));
}

/**
 * Actual/Actual (ISDA), §4.16(b). Split the period at each 1-January and divide each piece by the
 * length of *its own* year:
 *
 *   τ = Σ_y  days_in_period_within_year_y / (365 or 366 for year y)
 *
 * Equivalently, for start in year Y₁ and end in year Y₂ > Y₁:
 *   τ = (Jan1(Y₁+1) − start)/len(Y₁) + (Y₂ − Y₁ − 1) + (end − Jan1(Y₂))/len(Y₂)
 * and τ = (end − start)/len(Y₁) when Y₁ = Y₂. Handle end < start by returning −τ(end, start).
 *
 * Golden (ISDA / EMU conventions memo): 2003-11-01 → 2004-05-01 = 61/365 + 121/366 = 0.497724.
 * Summing τ over consecutive periods is exact (telescopes), which ACT/365F does not satisfy across a
 * leap year — that is a property test.
 */
double yfActActIsda(Date start, Date end) {
    if (start > end) return -yfActActIsda(end, start); // antisymmetric, like every other convention here
    const int y1 = start.year();
    const int y2 = end.year();
    auto yearLength = [](int y) { return Date::isLeap(y) ? 366.0 : 365.0; };
    if (y1 == y2) return static_cast<double>(end - start) / yearLength(y1); // double division: int/int would truncate to 0
    Date jan1Next = Date(y1 + 1, 1, 1);
    Date jan1End = Date(y2, 1, 1);

    const double firstPiece = static_cast<double>(jan1Next - start) / yearLength(y1);
    const double middlePiece = static_cast<double>(y2 - y1 - 1);
    const double lastPiece = static_cast<double>(end - jan1End) / yearLength(y2);

    return firstPiece + middlePiece + lastPiece;
}

} // namespace detail
} // namespace rl::core
