#pragma once
// Business-day calendars and conventions.
//
// A Calendar is a name plus a sorted list of holiday serials loaded from data/holidays/<file>.txt,
// plus the Saturday/Sunday weekend rule. It never computes a holiday: production systems take
// holidays from reference data, and tools/gen_holidays.py plays that role here. Calendars are
// immutable, cheap to copy (shared_ptr to the holiday table), and thread-safe to read.
//
// Registry names (Calendar::get): "USD-SIFMA", "USD-FEDWIRE", "GBP", "EUR-TARGET2", "JPY", "NONE".
// Joint calendars: "USD-SIFMA+GBP" (holiday if a holiday in either).

#include "rl/core/date.hpp"
#include "rl/core/tenor.hpp"

#include <filesystem>
#include <memory>
#include <string>
#include <string_view>
#include <vector>

namespace rl::core {

enum class BusinessDayConvention {
    Unadjusted,         ///< Leave the date alone (accrual dates on some bond schedules).
    Following,          ///< Next business day.
    ModifiedFollowing,  ///< Next business day unless that crosses a month end, then preceding. Swaps.
    Preceding,          ///< Previous business day.
    ModifiedPreceding,  ///< Previous business day unless that crosses a month start, then following.
};
std::string_view toString(BusinessDayConvention bdc);
/// "U"/"UNADJUSTED", "F"/"FOLLOWING", "MF"/"MODFOLLOWING"/"MODIFIEDFOLLOWING", "P", "MP". Throws.
BusinessDayConvention bdcFromString(std::string_view s);

/// Directory holding holidays/ and sobol/: env RL_DATA_DIR if set, else the configure-time path.
std::filesystem::path dataDir();

class Calendar {
public:
    /// Weekends-only calendar named "NONE".
    Calendar();

    /// Lookup by registry name (see file header). Loads and caches the holiday file on first use.
    static Calendar get(std::string_view name);
    static Calendar weekendsOnly() { return Calendar(); }
    /// Build from an explicit holiday list (tests, ad-hoc calendars).
    static Calendar fromHolidays(std::string name, std::vector<Date> holidays);
    /// Parse a holiday file: one "YYYY-MM-DD" per line, '#' comments, blank lines ignored.
    static Calendar fromFile(std::string name, const std::filesystem::path& file);
    /// Union of holidays: a day is a holiday if it is a holiday in either. Name "A+B".
    static Calendar join(const Calendar& a, const Calendar& b);

    const std::string& name() const noexcept { return impl_->name; }
    const std::vector<Date>& holidays() const noexcept { return impl_->holidays; }
    /// First/last holiday in the table; dates outside this range are weekends-only (loud in README).
    Date firstHoliday() const;
    Date lastHoliday() const;

    bool isWeekend(Date d) const { return d.isWeekend(); }
    bool isHoliday(Date d) const; // binary search on the sorted table; weekends are not "holidays"
    bool isBusinessDay(Date d) const { return !isWeekend(d) && !isHoliday(d); }

    // --- the maths of this module (calendar.cpp) ----------------------------------------------
    Date adjust(Date d, BusinessDayConvention bdc) const;
    /// Move by n business days (n < 0 moves back; n == 0 returns d unchanged even if non-business).
    Date advance(Date d, int businessDays) const;
    /// d + tenor, then adjust. For month/year tenors with eom=true and d the last business day of
    /// its month, the result is rolled to the last business day of the target month.
    Date advance(Date d, Tenor t, BusinessDayConvention bdc, bool eom = false) const;
    /// Number of business days in [from, to) by default; flags include/exclude the endpoints.
    int businessDaysBetween(Date from, Date to, bool includeFrom = true, bool includeTo = false) const;

    friend bool operator==(const Calendar& a, const Calendar& b) noexcept { return a.impl_ == b.impl_ || a.name() == b.name(); }

private:
    struct Impl {
        std::string name;
        std::vector<Date> holidays; // sorted, unique
    };
    explicit Calendar(std::shared_ptr<const Impl> impl) : impl_(std::move(impl)) {}
    std::shared_ptr<const Impl> impl_;
};

} // namespace rl::core
