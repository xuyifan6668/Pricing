#pragma once
// Swap: id + a fixed leg + a floating leg with opposite sides and one currency. Which one is paid is
// on the legs; payer/receiver is derived. Value type; owns its legs.
//
// DESIGN: no common Instrument base class. Per-trade dispatch  is a registry on a type tag /
// std::variant of trade types; a virtual hierarchy here would only exist to be
// downcast from.

#include "rl/core/date.hpp"
#include "rl/instruments/leg.hpp"

#include <string>

namespace rl::instruments {

class Swap {
public:
    /// Throws InputError unless the legs have opposite sides and the same currency.
    Swap(std::string id, Leg fixed, Leg floating);

    const std::string& id() const noexcept { return id_; }
    const Leg& fixedLeg() const noexcept { return fixed_; }
    const Leg& floatLeg() const noexcept { return float_; }
    bool payerOfFixed() const noexcept { return fixed_.side() == Side::Pay; }
    core::Currency currency() const noexcept { return fixed_.currency(); }
    core::Date maturity() const; ///< the later of the two legs' last payments

private:
    std::string id_;
    Leg fixed_;
    Leg float_;
};

} // namespace rl::instruments
