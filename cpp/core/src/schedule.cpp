#include "rl/core/schedule.hpp"

#include "rl/core/error.hpp"

#include <algorithm>
#include <cctype>
#include <string>
#include <utility>

namespace rl::core {

// ------------------------------------------------------------------------------------------------
// Plumbing
// ------------------------------------------------------------------------------------------------

std::string_view toString(StubRule s) {
    switch (s) {
    case StubRule::None: return "None";
    case StubRule::ShortFront: return "ShortFront";
    case StubRule::LongFront: return "LongFront";
    case StubRule::ShortBack: return "ShortBack";
    case StubRule::LongBack: return "LongBack";
    }
    return "?";
}

StubRule stubRuleFromString(std::string_view in) {
    std::string s;
    for (char c : in)
        if (std::isalnum(static_cast<unsigned char>(c))) s.push_back(static_cast<char>(std::toupper(static_cast<unsigned char>(c))));
    if (s == "NONE") return StubRule::None;
    if (s == "SHORTFRONT" || s == "SF") return StubRule::ShortFront;
    if (s == "LONGFRONT" || s == "LF") return StubRule::LongFront;
    if (s == "SHORTBACK" || s == "SB") return StubRule::ShortBack;
    if (s == "LONGBACK" || s == "LB") return StubRule::LongBack;
    throw InputError("unknown stub rule: " + std::string(in));
}

Schedule Schedule::generate(const ScheduleSpec& spec) { return build(spec, generateUnadjusted(spec)); }

Schedule Schedule::build(const ScheduleSpec& spec, std::vector<Date> unadj) {
    RL_REQUIRE(unadj.size() >= 2, "schedule needs at least two dates");
    RL_REQUIRE(unadj.front() == spec.effective, "first unadjusted date must be the effective date");
    RL_REQUIRE(unadj.back() == spec.termination, "last unadjusted date must be the termination date");
    for (std::size_t i = 1; i < unadj.size(); ++i)
        RL_REQUIRE(unadj[i] > unadj[i - 1], "unadjusted schedule dates must be strictly increasing");

    const Calendar& payCal = spec.paymentCalendar ? *spec.paymentCalendar : spec.calendar;
    const std::size_t n = unadj.size() - 1;

    // A period is "regular" if its length is exactly the frequency measured from its own start
    // (month arithmetic with clamping, or calendar days for D/W tenors). Only the first and last
    // period can be irregular; everything between is regular by construction.
    auto isRegular = [&](Date s, Date e) {
        if (spec.frequency.isMonthly()) {
            // Both the exact roll and the EOM-clamped roll count as regular.
            Date exact = s.addMonths(spec.frequency.inMonths());
            return e == exact || (spec.endOfMonth && e == exact.endOfMonth());
        }
        return e - s == spec.frequency.inDays();
    };

    std::vector<Period> periods;
    periods.reserve(n);
    for (std::size_t i = 0; i < n; ++i) {
        Period p;
        p.unadjustedStart = unadj[i];
        p.unadjustedEnd = unadj[i + 1];
        p.start = spec.calendar.adjust(p.unadjustedStart, spec.convention);
        p.end = spec.calendar.adjust(p.unadjustedEnd, i + 1 == n ? spec.terminationConvention : spec.convention);
        p.payment = spec.paymentLag == 0 ? p.end : payCal.advance(p.end, spec.paymentLag);
        p.regular = (i == 0 || i + 1 == n) ? isRegular(p.unadjustedStart, p.unadjustedEnd) : true;
        periods.push_back(p);
    }
    return Schedule(spec, std::move(periods));
}

std::vector<Date> Schedule::dates() const {
    std::vector<Date> out;
    out.reserve(periods_.size() + 1);
    if (periods_.empty()) return out;
    out.push_back(periods_.front().start);
    for (const auto& p : periods_) out.push_back(p.end);
    return out;
}

std::vector<Date> Schedule::unadjustedDates() const {
    std::vector<Date> out;
    out.reserve(periods_.size() + 1);
    if (periods_.empty()) return out;
    out.push_back(periods_.front().unadjustedStart);
    for (const auto& p : periods_) out.push_back(p.unadjustedEnd);
    return out;
}

std::vector<Date> Schedule::paymentDates() const {
    std::vector<Date> out;
    out.reserve(periods_.size());
    for (const auto& p : periods_) out.push_back(p.payment);
    return out;
}

// ------------------------------------------------------------------------------------------------
// Unadjusted date generation
// ------------------------------------------------------------------------------------------------

/**
 * Unadjusted period boundaries (ISDA 2006 §4.10–4.12 vocabulary; see the header).
 *
 * Generation is always from an anchor with anchor + k·frequency (never iteratively from the previous
 * clamped date: 31-Jan → 28-Feb → 28-Mar would lose the 31). The anchor and direction follow the stub
 * rule: front stubs and None generate backwards from termination (or nextToLastDate); back stubs
 * generate forwards from effective (or firstRegularDate).
 *
 * Roll / EOM, applied to every generated date but never to effective, termination or the explicit
 * dates: with endOfMonth and an end-of-month anchor every generated date is its month end; else with
 * rollDay the day is min(rollDay, daysInMonth); else whatever addMonths gives. Day/week frequencies
 * add calendar days and ignore roll/EOM.
 *
 * A stub exists iff the generated sequence does not land exactly on the far boundary. LongFront /
 * LongBack merge the stub into the neighbouring regular period only when a stub exists; None throws
 * when one does.
 */
std::vector<Date> Schedule::generateUnadjusted(const ScheduleSpec& spec) {
    RL_REQUIRE(spec.effective < spec.termination, "schedule: effective must precede termination");
    RL_REQUIRE(spec.frequency.n > 0, "schedule: frequency must be positive");
    if (spec.rollDay) RL_REQUIRE(*spec.rollDay >= 1 && *spec.rollDay <= 31, "schedule: rollDay must be in [1, 31]");
    if (spec.firstRegularDate)
        RL_REQUIRE(*spec.firstRegularDate > spec.effective && *spec.firstRegularDate < spec.termination,
                   "schedule: firstRegularDate must lie strictly inside (effective, termination)");
    if (spec.nextToLastDate)
        RL_REQUIRE(*spec.nextToLastDate > spec.effective && *spec.nextToLastDate < spec.termination,
                   "schedule: nextToLastDate must lie strictly inside (effective, termination)");

    const bool monthly = spec.frequency.isMonthly();

    // anchor + k periods, k signed.
    auto hop = [&](Date anchor, int k) {
        return monthly ? anchor.addMonths(spec.frequency.inMonths() * k) : anchor.addDays(spec.frequency.inDays() * k);
    };
    // Roll-day / end-of-month rule on a generated date.
    auto roll = [&](Date raw, Date anchor) {
        if (!monthly) return raw;
        if (spec.endOfMonth && anchor.isEndOfMonth()) return raw.endOfMonth();
        if (spec.rollDay) return Date(raw.year(), raw.month(), std::min(static_cast<unsigned>(*spec.rollDay), Date::daysInMonth(raw.year(), raw.month())));
        return raw;
    };

    // Walk from `anchor` towards `boundary` (sign = −1 backwards to effective, +1 forwards to termination).
    // Returns the generated dates strictly inside (effective, termination) in walking order, and
    // whether the walk landed exactly on the boundary (no stub).
    auto walk = [&](Date anchor, int sign, Date boundary) {
        std::vector<Date> out;
        bool exact = false;
        for (int k = 1;; ++k) {
            const Date d = roll(hop(anchor, sign * k), anchor);
            if (d == boundary) {
                exact = true;
                break;
            }
            if (sign < 0 ? d < boundary : d > boundary) break;
            out.push_back(d);
        }
        return std::pair{out, exact};
    };

    std::vector<Date> dates{spec.effective};
    const bool frontRule = spec.stub == StubRule::ShortFront || spec.stub == StubRule::LongFront || spec.stub == StubRule::None;

    if (spec.firstRegularDate) {
        // Explicit front stub; regular periods run forward from firstRegularDate.
        dates.push_back(*spec.firstRegularDate);
        auto [gen, exact] = walk(*spec.firstRegularDate, +1, spec.nextToLastDate.value_or(spec.termination));
        if (spec.stub == StubRule::LongBack && !exact && !gen.empty()) gen.pop_back();
        dates.insert(dates.end(), gen.begin(), gen.end());
    } else if (frontRule) {
        const Date anchor = spec.nextToLastDate.value_or(spec.termination);
        auto [gen, exact] = walk(anchor, -1, spec.effective); // descending
        if (spec.stub == StubRule::None) RL_REQUIRE(exact, "schedule: StubRule::None but the periods do not land on the effective date");
        if (spec.stub == StubRule::LongFront && !exact && !gen.empty()) gen.pop_back(); // merge the short stub into the first regular period
        dates.insert(dates.end(), gen.rbegin(), gen.rend());
    } else {
        auto [gen, exact] = walk(spec.effective, +1, spec.nextToLastDate.value_or(spec.termination));
        if (spec.stub == StubRule::LongBack && !exact && !gen.empty()) gen.pop_back();
        dates.insert(dates.end(), gen.begin(), gen.end());
    }
    if (spec.nextToLastDate && dates.back() != *spec.nextToLastDate) dates.push_back(*spec.nextToLastDate);
    dates.push_back(spec.termination);
    return dates;
}

} // namespace rl::core
