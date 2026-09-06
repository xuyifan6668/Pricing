#include "rl/core/tenor.hpp"

#include "rl/core/error.hpp"

#include <cctype>
#include <charconv>
#include <ostream>

namespace rl::core {

Tenor Tenor::parse(std::string_view in) {
    std::string s;
    for (char c : in)
        if (!std::isspace(static_cast<unsigned char>(c))) s.push_back(static_cast<char>(std::toupper(static_cast<unsigned char>(c))));
    RL_REQUIRE(!s.empty(), "empty tenor string");
    if (s == "ON" || s == "O/N") return days(1);
    if (s == "TN" || s == "T/N") return days(2);
    if (s == "SN" || s == "S/N") return days(3);
    char u = s.back();
    int n = 0;
    auto [p, ec] = std::from_chars(s.data(), s.data() + s.size() - 1, n);
    RL_REQUIRE(ec == std::errc{} && p == s.data() + s.size() - 1, "cannot parse tenor: " + std::string(in));
    switch (u) {
    case 'D': return days(n);
    case 'W': return weeks(n);
    case 'M': return months(n);
    case 'Y': return years(n);
    default: throw InputError("unknown tenor unit in: " + std::string(in));
    }
}

std::string Tenor::str() const {
    static constexpr char suffix[] = {'D', 'W', 'M', 'Y'};
    return std::to_string(n) + suffix[static_cast<int>(unit)];
}

int Tenor::inMonths() const {
    switch (unit) {
    case TenorUnit::Months: return n;
    case TenorUnit::Years: return 12 * n;
    default: throw InputError("tenor " + str() + " is not a month/year tenor");
    }
}

int Tenor::inDays() const {
    switch (unit) {
    case TenorUnit::Days: return n;
    case TenorUnit::Weeks: return 7 * n;
    default: throw InputError("tenor " + str() + " is not a day/week tenor");
    }
}

Tenor Tenor::normalised() const noexcept {
    if (unit == TenorUnit::Months && n % 12 == 0 && n != 0) return years(n / 12);
    if (unit == TenorUnit::Days && n % 7 == 0 && n != 0) return weeks(n / 7);
    return *this;
}

std::ostream& operator<<(std::ostream& os, Tenor t) { return os << t.str(); }

std::string_view toString(TenorUnit u) {
    static constexpr std::string_view names[] = {"Days", "Weeks", "Months", "Years"};
    return names[static_cast<int>(u)];
}

} // namespace rl::core
