#pragma once
// Serial-day date. One int, no time zone, no time of day: a rates library never needs either.
//
// Epoch [D11]: serial 0 == 1899-12-30, so serial() equals the Excel/OLE date number for every date
// on or after 1900-03-01 (Excel's fictitious 1900-02-29 is why the epoch is the 30th and not the
// 31st). A trader pasting a date from Excel into the add-in gets the same integer back, with no
// conversion layer to get wrong. 1970-01-01 is serial 25569.
//
// Everything here is passive (never differentiated) and trivially copyable; pass by value.

#include <compare>
#include <cstdint>
#include <iosfwd>
#include <string>
#include <string_view>

namespace rl::core {

enum class Weekday : int { Sunday = 0, Monday, Tuesday, Wednesday, Thursday, Friday, Saturday };
enum class Month : unsigned {
    January = 1, February, March, April, May, June, July, August, September, October, November, December
};

struct Ymd {
    int year;
    unsigned month; // 1..12
    unsigned day;   // 1..31
    friend bool operator==(const Ymd&, const Ymd&) noexcept = default;
};

class Date {
public:
    /// Days between the Unix epoch (1970-01-01) and our epoch (1899-12-30).
    static constexpr int kUnixEpochSerial = 25569;

    constexpr Date() noexcept : serial_(0) {}
    constexpr explicit Date(int serial) noexcept : serial_(serial) {}
    /// Validates the civil triple (throws InputError) then converts.
    Date(int year, unsigned month, unsigned day);
    Date(int year, Month month, unsigned day) : Date(year, static_cast<unsigned>(month), day) {}

    static constexpr Date fromSerial(int serial) noexcept { return Date(serial); }
    static Date fromYmd(int year, unsigned month, unsigned day) { return Date(year, month, day); }
    /// Parses "YYYY-MM-DD" (ISO 8601 extended); throws InputError on anything else.
    static Date fromIso(std::string_view iso);
    static constexpr Date minDate() noexcept { return Date(-693593); } // 0001-01-01
    static constexpr Date maxDate() noexcept { return Date(2958465); } // 9999-12-31

    constexpr int serial() const noexcept { return serial_; }
    Ymd ymd() const { return civilFromDays(serial_); }
    int year() const { return ymd().year; }
    unsigned month() const { return ymd().month; }
    unsigned day() const { return ymd().day; }
    Weekday weekday() const;
    bool isEndOfMonth() const { return day() == daysInMonth(year(), month()); }
    bool isWeekend() const {
        auto w = weekday();
        return w == Weekday::Saturday || w == Weekday::Sunday;
    }
    Date endOfMonth() const;
    Date addDays(int n) const noexcept { return Date(serial_ + n); }
    Date addMonths(int n) const;
    Date addYears(int n) const { return addMonths(12 * n); }
    /// "YYYY-MM-DD"
    std::string toIso() const;

    // --- calendar arithmetic (the maths of this module; see date.cpp) --------------------------
    static bool isLeap(int year);
    static unsigned daysInMonth(int year, unsigned month);
    /// Civil (proleptic Gregorian) -> serial. No validation; caller guarantees a valid triple.
    static int daysFromCivil(int year, unsigned month, unsigned day);
    /// Serial -> civil. Inverse of daysFromCivil for every int in [minDate, maxDate].
    static Ymd civilFromDays(int serial);

    // --- operators ----------------------------------------------------------------------------
    // DESIGN: a defaulted <=> also synthesises == (a user-provided <=> does not), so this one line
    // gives all six comparisons; a hand-written <=> silently drops == and every `d1 == d2` stops compiling.
    friend constexpr auto operator<=>(Date a, Date b) noexcept = default;
    friend constexpr Date operator+(Date d, int days) noexcept { return Date(d.serial_ + days); }
    friend constexpr Date operator-(Date d, int days) noexcept { return Date(d.serial_ - days); }
    /// Actual days from b to a (a - b), signed.
    friend constexpr int operator-(Date a, Date b) noexcept { return a.serial_ - b.serial_; }
    Date& operator+=(int days) noexcept { serial_ += days; return *this; }
    Date& operator-=(int days) noexcept { serial_ -= days; return *this; }
    Date& operator++() noexcept { ++serial_; return *this; }
    Date& operator--() noexcept { --serial_; return *this; }

private:
    int serial_;
};

std::ostream& operator<<(std::ostream& os, Date d);
std::string_view toString(Weekday w);
std::string_view toString(Month m);

} // namespace rl::core
