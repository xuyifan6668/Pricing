#pragma once
// Tenor: an integer count of days, weeks, months or years. Pure description; the calendar decides
// what "3M from today" means (see Calendar::advance).

#include <compare>
#include <iosfwd>
#include <string>
#include <string_view>

namespace rl::core {

enum class TenorUnit { Days, Weeks, Months, Years };

struct Tenor {
    int n = 0;
    TenorUnit unit = TenorUnit::Days;

    static constexpr Tenor days(int n) noexcept { return {n, TenorUnit::Days}; }
    static constexpr Tenor weeks(int n) noexcept { return {n, TenorUnit::Weeks}; }
    static constexpr Tenor months(int n) noexcept { return {n, TenorUnit::Months}; }
    static constexpr Tenor years(int n) noexcept { return {n, TenorUnit::Years}; }

    /// "3M", "10Y", "1W", "7D", "ON" (=1D), "TN" (=2D), "SN"(=3D); case-insensitive; throws InputError.
    static Tenor parse(std::string_view s);
    /// Canonical form: "3M", "10Y", "1W", "7D".
    std::string str() const;

    /// True if the unit is months or years (i.e. month-arithmetic with day clamping applies).
    constexpr bool isMonthly() const noexcept { return unit == TenorUnit::Months || unit == TenorUnit::Years; }
    /// Months for month/year tenors; throws InputError for day/week tenors.
    int inMonths() const;
    /// Calendar days for day/week tenors; throws InputError for month/year tenors.
    int inDays() const;
    /// Normalises 12M -> 1Y, 7D -> 1W; used for equality of "the same" tenor written differently.
    Tenor normalised() const noexcept;

    friend constexpr Tenor operator*(Tenor t, int k) noexcept { return {t.n * k, t.unit}; }
    friend constexpr Tenor operator*(int k, Tenor t) noexcept { return {t.n * k, t.unit}; }
    friend constexpr Tenor operator-(Tenor t) noexcept { return {-t.n, t.unit}; }
    friend bool operator==(Tenor a, Tenor b) noexcept { return a.normalised().n == b.normalised().n && a.normalised().unit == b.normalised().unit; }
};

std::ostream& operator<<(std::ostream& os, Tenor t);
std::string_view toString(TenorUnit u);

} // namespace rl::core
