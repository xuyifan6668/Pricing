#include "rl/core/calendar.hpp"

#include "rl/core/error.hpp"

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <fstream>
#include <map>
#include <mutex>

#ifndef RL_DATA_DIR
#define RL_DATA_DIR "data"
#endif

namespace rl::core {

// ------------------------------------------------------------------------------------------------
// Plumbing
// ------------------------------------------------------------------------------------------------

std::string_view toString(BusinessDayConvention bdc) {
    switch (bdc) {
    case BusinessDayConvention::Unadjusted: return "Unadjusted";
    case BusinessDayConvention::Following: return "Following";
    case BusinessDayConvention::ModifiedFollowing: return "ModifiedFollowing";
    case BusinessDayConvention::Preceding: return "Preceding";
    case BusinessDayConvention::ModifiedPreceding: return "ModifiedPreceding";
    }
    return "?";
}

BusinessDayConvention bdcFromString(std::string_view in) {
    std::string s;
    for (char c : in)
        if (std::isalnum(static_cast<unsigned char>(c))) s.push_back(static_cast<char>(std::toupper(static_cast<unsigned char>(c))));
    if (s == "U" || s == "UNADJUSTED" || s == "NONE") return BusinessDayConvention::Unadjusted;
    if (s == "F" || s == "FOLLOWING") return BusinessDayConvention::Following;
    if (s == "MF" || s == "MODFOLLOWING" || s == "MODIFIEDFOLLOWING") return BusinessDayConvention::ModifiedFollowing;
    if (s == "P" || s == "PRECEDING") return BusinessDayConvention::Preceding;
    if (s == "MP" || s == "MODPRECEDING" || s == "MODIFIEDPRECEDING") return BusinessDayConvention::ModifiedPreceding;
    throw InputError("unknown business-day convention: " + std::string(in));
}

std::filesystem::path dataDir() {
    if (const char* env = std::getenv("RL_DATA_DIR"); env && *env) return env;
    return RL_DATA_DIR;
}

Calendar::Calendar() : impl_(std::make_shared<const Impl>(Impl{"NONE", {}})) {}

Calendar Calendar::fromHolidays(std::string name, std::vector<Date> holidays) {
    std::sort(holidays.begin(), holidays.end());
    holidays.erase(std::unique(holidays.begin(), holidays.end()), holidays.end());
    // DESIGN: weekend entries are kept as loaded. holidays() is the reference-data table, and a
    // holiday file that lists 1 January on a Saturday must round-trip unchanged (the loader is not a
    // second source of truth for what a closure is; isBusinessDay() combines the two rules).
    return Calendar(std::make_shared<const Impl>(Impl{std::move(name), std::move(holidays)}));
}

Calendar Calendar::fromFile(std::string name, const std::filesystem::path& file) {
    std::ifstream in(file);
    if (!in) throw InputError("cannot open holiday file: " + file.string());
    std::vector<Date> hols;
    std::string line;
    int lineNo = 0;
    while (std::getline(in, line)) {
        ++lineNo;
        auto hash = line.find('#');
        if (hash != std::string::npos) line.erase(hash);
        auto b = line.find_first_not_of(" \t\r");
        if (b == std::string::npos) continue;
        auto e = line.find_last_not_of(" \t\r");
        std::string tok = line.substr(b, e - b + 1);
        hols.push_back(withContext(file.string() + ":" + std::to_string(lineNo), [&] { return Date::fromIso(tok); }));
    }
    return fromHolidays(std::move(name), std::move(hols));
}

namespace {
const std::map<std::string, std::string, std::less<>>& registryFiles() {
    static const std::map<std::string, std::string, std::less<>> files = {
        {"USD-SIFMA", "us_sifma.txt"}, {"USD-FEDWIRE", "us_fedwire.txt"}, {"GBP", "uk.txt"},
        {"EUR-TARGET2", "target2.txt"}, {"JPY", "tokyo.txt"},
    };
    return files;
}
} // namespace

Calendar Calendar::get(std::string_view name) {
    if (name == "NONE" || name.empty()) return Calendar();
    if (auto plus = name.find('+'); plus != std::string_view::npos)
        return join(get(name.substr(0, plus)), get(name.substr(plus + 1)));

    static std::mutex mtx;
    static std::map<std::string, Calendar, std::less<>> cache;
    std::lock_guard lock(mtx);
    if (auto it = cache.find(name); it != cache.end()) return it->second;

    auto f = registryFiles().find(name);
    if (f == registryFiles().end()) throw InputError("unknown calendar: " + std::string(name));
    Calendar cal = fromFile(std::string(name), dataDir() / "holidays" / f->second);
    cache.emplace(std::string(name), cal);
    return cal;
}

Calendar Calendar::join(const Calendar& a, const Calendar& b) {
    std::vector<Date> hols(a.holidays());
    hols.insert(hols.end(), b.holidays().begin(), b.holidays().end());
    return fromHolidays(a.name() + "+" + b.name(), std::move(hols));
}

Date Calendar::firstHoliday() const {
    RL_REQUIRE(!impl_->holidays.empty(), "calendar " + name() + " has no holidays");
    return impl_->holidays.front();
}

Date Calendar::lastHoliday() const {
    RL_REQUIRE(!impl_->holidays.empty(), "calendar " + name() + " has no holidays");
    return impl_->holidays.back();
}

bool Calendar::isHoliday(Date d) const {
    const auto& h = impl_->holidays;
    return std::binary_search(h.begin(), h.end(), d);
}

// ------------------------------------------------------------------------------------------------
// Business-day arithmetic
// ------------------------------------------------------------------------------------------------

namespace {
// DESIGN: adjustment loops are unbounded on purpose. The longest closure run in any v1 calendar is a
// weekend plus a holiday cluster (Tokyo Golden Week, ~6 days); a loop capped at N days would return a
// non-business day silently once a calendar with a longer run appears, which is worse than looping.
Date following(const Calendar& cal, Date d) {
    while (!cal.isBusinessDay(d)) d = d.addDays(1);
    return d;
}
Date preceding(const Calendar& cal, Date d) {
    while (!cal.isBusinessDay(d)) d = d.addDays(-1);
    return d;
}
} // namespace

/**
 * Apply a business-day convention (ISDA 2006 §4.12):
 *   Unadjusted:         d
 *   Following:          smallest business day ≥ d
 *   Preceding:          largest business day ≤ d
 *   ModifiedFollowing:  Following, unless that lands in a later month, then Preceding
 *   ModifiedPreceding:  Preceding, unless that lands in an earlier month, then Following
 * A business day is a fixed point under every convention.
 */
Date Calendar::adjust(Date d, BusinessDayConvention bdc) const {
    switch (bdc) {
    case BusinessDayConvention::Unadjusted: return d;
    case BusinessDayConvention::Following: return following(*this, d);
    case BusinessDayConvention::Preceding: return preceding(*this, d);
    case BusinessDayConvention::ModifiedFollowing: {
        const Date f = following(*this, d);
        return f.month() == d.month() ? f : preceding(*this, d);
    }
    case BusinessDayConvention::ModifiedPreceding: {
        const Date p = preceding(*this, d);
        return p.month() == d.month() ? p : following(*this, d);
    }
    }
    throw InputError("Calendar::adjust: unknown BusinessDayConvention");
}

/**
 * Step n business days: n > 0 moves forward n business days, n < 0 backward; n == 0 returns d
 * unchanged even if d is not a business day (callers wanting "today or next" use adjust(Following)).
 * This is the settlement-lag operation: spot = advance(today, 2); SOFR payment = advance(end, 2).
 */
Date Calendar::advance(Date d, int businessDays) const {
    const int step = businessDays > 0 ? 1 : -1;
    for (int i = 0; i < std::abs(businessDays); ++i) {
        d = d.addDays(step);
        while (!isBusinessDay(d)) d = d.addDays(step);
    }
    return d;
}

/**
 * d + tenor, then adjust(bdc). Day/week tenors add calendar days (*not* business days — a 1W deposit
 * is 7 calendar days then adjusted). Month/year tenors use Date::addMonths (day clamped).
 *
 * End-of-month rule (eom = true, month/year tenors only): if d is the last *business* day of its
 * month, the unadjusted result is the last calendar day of the target month (then adjusted). This is
 * how a 6M swap from 2026-02-27 (last business day of Feb) ends 2026-08-31, not 2026-08-27.
 * Test: advance(2026-02-27, 6M, MF, eom=true) under USD-SIFMA = 2026-08-31; with eom=false = 2026-08-27.
 */
Date Calendar::advance(Date d, Tenor t, BusinessDayConvention bdc, bool eom) const {
    Date unadjusted;
    switch (t.unit) {
    case TenorUnit::Days:
        unadjusted = d.addDays(t.n);
        break;
    case TenorUnit::Weeks:
        unadjusted = d.addDays(t.n * 7);
        break;
    case TenorUnit::Months:
        unadjusted = d.addMonths(t.n);
        break;
    case TenorUnit::Years:
        unadjusted = d.addMonths(t.n * 12);
        break;
    }
    // EOM rule: only for month/year tenors, and only when d is the last *business* day of its month.
    if (eom && t.isMonthly() && d == adjust(d.endOfMonth(), BusinessDayConvention::Preceding)) unadjusted = unadjusted.endOfMonth();
    return adjust(unadjusted, bdc);
}

/**
 * Count business days between two dates. Default counts the half-open [from, to): that is the
 * number of overnight fixings that accrue from `from` to `to`, the quantity an OIS coupon needs.
 * Return a negative count if to < from (swap and negate). O(to − from); fine for coupon periods.
 * Test: [2026-09-04, 2026-09-08) = 1 (Fri 4 only; Sat, Sun, Labor Day excluded).
 */
int Calendar::businessDaysBetween(Date from, Date to, bool includeFrom, bool includeTo) const {
    if (from > to) return -businessDaysBetween(to, from, includeTo, includeFrom);
    if (from == to) return (includeFrom && includeTo && isBusinessDay(from)) ? 1 : 0;
    int count = 0;
    if (includeFrom && isBusinessDay(from)) ++count;
    if (includeTo && isBusinessDay(to)) ++count;
    for (Date d = from + 1; d < to; ++d)
        if (isBusinessDay(d)) ++count;
    return count;
}

} // namespace rl::core
