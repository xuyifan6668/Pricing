#pragma once
// BasisSwap: two floating legs (overnight vs overnight in v1: FF + spread vs SOFR flat), opposite
// sides, one currency. legs()[0] is the spread leg by builder convention; the spread itself is on the
// coupons of that leg (all equal), not on the trade, so the leg PV code has nothing special to know.

#include "rl/core/date.hpp"
#include "rl/instruments/leg.hpp"

#include <array>
#include <string>

namespace rl::instruments {

class BasisSwap {
public:
    /// Throws InputError unless the legs have opposite sides and the same currency.
    BasisSwap(std::string id, Leg first, Leg second);

    const std::string& id() const noexcept { return id_; }
    const std::array<Leg, 2>& legs() const noexcept { return legs_; }
    core::Currency currency() const noexcept { return legs_[0].currency(); }
    core::Date maturity() const;

private:
    std::string id_;
    std::array<Leg, 2> legs_;
};

} // namespace rl::instruments
