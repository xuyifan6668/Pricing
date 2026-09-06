#pragma once
// Currencies in v1 scope (USD is the only one priced; the rest exist so Index and conventions
// serialise correctly and v2 is a data change, not a code change).

#include <string_view>

namespace rl::core {

enum class Currency { USD, EUR, GBP, JPY };

std::string_view toString(Currency c);
Currency currencyFromString(std::string_view s); // throws InputError

/// Spot lag in business days for cash-settled instruments in that currency.
constexpr int settlementDays(Currency c) noexcept {
    switch (c) {
    case Currency::USD: return 2;
    case Currency::EUR: return 2;
    case Currency::GBP: return 0;
    case Currency::JPY: return 2;
    }
    return 2;
}

/// Registry name of the settlement calendar (see Calendar::get).
constexpr std::string_view settlementCalendarName(Currency c) noexcept {
    switch (c) {
    case Currency::USD: return "USD-FEDWIRE";
    case Currency::EUR: return "EUR-TARGET2";
    case Currency::GBP: return "GBP";
    case Currency::JPY: return "JPY";
    }
    return "NONE";
}

} // namespace rl::core
