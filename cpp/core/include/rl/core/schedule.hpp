#pragma once
// Coupon schedule generation with stubs, end-of-month and roll-day handling.
//
// Vocabulary (ISDA 2006 §4.10 – §4.12, and every confirm you will ever read):
//   effective date       first accrual start (unadjusted; usually spot)
//   termination date     last accrual end (unadjusted)
//   frequency            3M / 6M / 1Y between regular period dates
//   stub                 the one irregular period, at the front or the back, shorter or longer than
//                        a regular period. ShortFront is the market default for swaps.
//   roll day             the day-of-month every regular date falls on (before adjustment). Comes
//                        from the anchor date (termination for front stubs, effective for back
//                        stubs) unless given explicitly.
//   EOM rule             if the anchor is the last day of its month, every regular date is the last
//                        day of its month (so 30-Nov rolls to 28/29-Feb, 31-May, 31-Aug).
//   payment lag          business days from adjusted accrual end to payment (SOFR OIS: 2).

#include "rl/core/calendar.hpp"
#include "rl/core/date.hpp"
#include "rl/core/tenor.hpp"

#include <cstddef>
#include <optional>
#include <string_view>
#include <vector>

namespace rl::core {

enum class StubRule { None, ShortFront, LongFront, ShortBack, LongBack };
std::string_view toString(StubRule s);
StubRule stubRuleFromString(std::string_view s);

struct ScheduleSpec {
    Date effective;
    Date termination;
    Tenor frequency = Tenor::months(6);
    Calendar calendar = Calendar::weekendsOnly();
    BusinessDayConvention convention = BusinessDayConvention::ModifiedFollowing;
    BusinessDayConvention terminationConvention = BusinessDayConvention::ModifiedFollowing;
    StubRule stub = StubRule::ShortFront;
    bool endOfMonth = false;
    /// Force the regular dates onto this day of month (1..31; 31 means "last day"). Overrides the
    /// anchor's day. Rarely needed once endOfMonth is set correctly; IMM-style rolls are v2.
    std::optional<int> rollDay;
    /// Explicit first regular period start (front stub end). Overrides `stub` for the front.
    std::optional<Date> firstRegularDate;
    /// Explicit last regular period end (back stub start). Overrides `stub` for the back.
    std::optional<Date> nextToLastDate;
    int paymentLag = 0;
    std::optional<Calendar> paymentCalendar; ///< Defaults to `calendar`.
};

struct Period {
    Date unadjustedStart;
    Date unadjustedEnd;
    Date start;   ///< adjusted accrual start
    Date end;     ///< adjusted accrual end
    Date payment; ///< end + paymentLag business days on the payment calendar
    bool regular; ///< false for the stub period(s)
};

class Schedule {
public:
    /// The whole algorithm. Equivalent to build(spec, generateUnadjusted(spec)).
    static Schedule generate(const ScheduleSpec& spec);

    /// Unadjusted period boundary dates, effective first, termination last, strictly increasing.
    /// The maths of this module (schedule.cpp).
    static std::vector<Date> generateUnadjusted(const ScheduleSpec& spec);

    /// Plumbing: adjust every date, mark the stub period(s), compute payment dates.
    /// `unadjusted` must be strictly increasing with front() == spec.effective, back() == spec.termination.
    static Schedule build(const ScheduleSpec& spec, std::vector<Date> unadjusted);

    const ScheduleSpec& spec() const noexcept { return spec_; }
    std::size_t size() const noexcept { return periods_.size(); }
    bool empty() const noexcept { return periods_.empty(); }
    const Period& operator[](std::size_t i) const noexcept { return periods_[i]; }
    const Period& front() const noexcept { return periods_.front(); }
    const Period& back() const noexcept { return periods_.back(); }
    const std::vector<Period>& periods() const noexcept { return periods_; }
    std::vector<Period>::const_iterator begin() const noexcept { return periods_.begin(); }
    std::vector<Period>::const_iterator end() const noexcept { return periods_.end(); }

    /// Adjusted boundary dates: size()+1 entries.
    std::vector<Date> dates() const;
    std::vector<Date> unadjustedDates() const;
    std::vector<Date> paymentDates() const;

private:
    Schedule(ScheduleSpec spec, std::vector<Period> periods) : spec_(std::move(spec)), periods_(std::move(periods)) {}
    ScheduleSpec spec_;
    std::vector<Period> periods_;
};

} // namespace rl::core
