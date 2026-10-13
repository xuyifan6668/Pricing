#pragma once
// Call/put, and the payer/receiver mapping for swaptions: a payer swaption is a call on the forward swap rate,
// a receiver swaption a put, both on the annuity numeraire (Andersen & Piterbarg vol. 1 §5.10). The formulas in
// this module are undiscounted (per unit annuity); the annuity is the caller's.

#include "rl/core/error.hpp"

#include <string>
#include <string_view>

namespace rl::vol {

enum class OptionType { Call, Put };

/// +1 for a call, −1 for a put: payoff = max(ω (F − K), 0).
constexpr double optionSign(OptionType t) noexcept { return t == OptionType::Call ? 1.0 : -1.0; }

constexpr OptionType optionTypeForPayer(bool payer) noexcept { return payer ? OptionType::Call : OptionType::Put; }

inline std::string_view toString(OptionType t) noexcept { return t == OptionType::Call ? "call" : "put"; }

/// "call" | "payer" → Call, "put" | "receiver" → Put. Throws InputError otherwise.
inline OptionType optionTypeFromString(std::string_view s) {
    if (s == "call" || s == "payer" || s == "c") return OptionType::Call;
    if (s == "put" || s == "receiver" || s == "p") return OptionType::Put;
    throw InputError("unknown option type '" + std::string(s) + "' (call|put|payer|receiver)");
}

} // namespace rl::vol
