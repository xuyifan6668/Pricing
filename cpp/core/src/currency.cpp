#include "rl/core/currency.hpp"

#include "rl/core/error.hpp"

#include <cctype>
#include <string>

namespace rl::core {

std::string_view toString(Currency c) {
    switch (c) {
    case Currency::USD: return "USD";
    case Currency::EUR: return "EUR";
    case Currency::GBP: return "GBP";
    case Currency::JPY: return "JPY";
    }
    return "?";
}

Currency currencyFromString(std::string_view in) {
    std::string s;
    for (char c : in) s.push_back(static_cast<char>(std::toupper(static_cast<unsigned char>(c))));
    if (s == "USD") return Currency::USD;
    if (s == "EUR") return Currency::EUR;
    if (s == "GBP") return Currency::GBP;
    if (s == "JPY") return Currency::JPY;
    throw InputError("unknown currency: " + std::string(in));
}

} // namespace rl::core
