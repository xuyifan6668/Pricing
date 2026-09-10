#pragma once
// Day-count conventions per ISDA 2006 Definitions §4.16. Enum + switch, not virtual (closed set,
// serialisable, fast inside the schedule builder — docs/DESIGN.md §2.3).
//
// Year fractions are passive doubles: a date difference is never differentiated.

#include "rl/core/date.hpp"

#include <string_view>

namespace rl::core {

enum class DayCount {
    Act360,      ///< ISDA 2006 §4.16(e) "Actual/360". SOFR, ESTR, Fed Funds, EURIBOR, USD swap floating legs.
    Act365Fixed, ///< §4.16(d) "Actual/365 (Fixed)". SONIA, TONA, GBP/JPY legs; also our curve-time convention.
    ActActIsda,  ///< §4.16(b) "Actual/Actual (ISDA)". Government bonds; the only one that splits at year ends.
    Thirty360US, ///< §4.16(f) "30/360" (Bond Basis). USD fixed swap legs.
    ThirtyE360,  ///< §4.16(g) "30E/360" (Eurobond Basis). EUR fixed swap legs.
};

/// "ACT/360", "ACT/365F", "ACT/ACT", "30/360", "30E/360"
std::string_view toString(DayCount dc);
/// Accepts the strings above plus common aliases ("A/360", "ACT/365", "ACT/ACT ISDA", "30/360 US",
/// "BOND", "30E/360 ISDA" is NOT accepted: it is a different convention). Throws InputError.
DayCount dayCountFromString(std::string_view s);

/// Day-count numerator: the (possibly convention-adjusted) number of days from start to end.
int dayCountDays(DayCount dc, Date start, Date end);
/// Accrual factor τ(start, end) under the convention. Negative if end < start (used by some
/// reference-period computations); callers that need non-negativity check it.
double yearFraction(DayCount dc, Date start, Date end);

namespace detail {
// One function per convention; the maths of this module. See daycount.cpp.
int days30360US(Date start, Date end);
int days30E360(Date start, Date end);
double yfActActIsda(Date start, Date end);
} // namespace detail

} // namespace rl::core
