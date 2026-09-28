#pragma once
// Single-period money-market instruments: Deposit, Fra, Future. Flat records rather than one-cashflow
// legs: their bootstrap residuals are rate identities, a future has no PV in the swap sense (daily
// margin), and forcing them into legs would invent a notional that does not exist on the contract.

#include "rl/core/date.hpp"
#include "rl/core/index.hpp"

#include <string>

namespace rl::instruments {

/// Cash deposit: lend `notional` from start to end at `rate` (simple, on the index's day count).
class Deposit {
public:
    Deposit(std::string id, core::Index index, core::Date start, core::Date end, double tau, double notional, double rate);

    const std::string& id() const noexcept { return id_; }
    const core::Index& index() const noexcept { return index_; }
    core::Date start() const noexcept { return start_; }
    core::Date end() const noexcept { return end_; }
    double tau() const noexcept { return tau_; }
    double notional() const noexcept { return notional_; }
    double rate() const noexcept { return rate_; }

private:
    std::string id_;
    core::Index index_;
    core::Date start_, end_;
    double tau_, notional_, rate_;
};

} // namespace rl::instruments
