#pragma once
// European swaption: id + the underlying swap (owned by value — a swaption is not a view on someone
// else's swap) + expiry + settlement type. Built and priced throughout; the shape is fixed now so
// that Phases 3–5 add pricers, not representation.
//
// Settlement: Physical (USD convention: exercise into the swap) vs Cash (EUR: pay the annuity-based
// cash amount); the difference matters for the annuity used and for CMS in v2.

#include "rl/core/date.hpp"
#include "rl/instruments/swap.hpp"

#include <string>

namespace rl::instruments {

class Swaption {
public:
    /// Throws InputError if expiry is after the underlying's first accrual start.
    Swaption(std::string id, Swap underlying, core::Date expiry, Settlement settlement);

    const std::string& id() const noexcept { return id_; }
    const Swap& underlying() const noexcept { return swap_; }
    core::Date expiry() const noexcept { return expiry_; }
    Settlement settlement() const noexcept { return settlement_; }
    /// Payer swaption = right to enter as the fixed-rate payer.
    bool payer() const noexcept { return swap_.payerOfFixed(); }

private:
    std::string id_;
    Swap swap_;
    core::Date expiry_;
    Settlement settlement_;
};

} // namespace rl::instruments
