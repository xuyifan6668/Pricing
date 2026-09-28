#pragma once
// Leg: an ordered, contiguous vector of cashflows with one side, one currency. Eager (materialised at
// build time), not lazy: the bootstrap iterates legs ~1e6 times per solve and a lazy leg would either
// regenerate its schedule or carry a cache (a mutex under the parallel risk driver). Memory is ~50 KB
// per 30y trade, which is nothing (docs/DESIGN.md "Trade representation").
//
// Invariant: cashflows are ordered by payment date (non-decreasing). The pay/receive sign lives here,
// not on the coupons: coupon notionals are positive.

#include "rl/core/currency.hpp"
#include "rl/core/date.hpp"
#include "rl/instruments/cashflow.hpp"

#include <vector>

namespace rl::instruments {

class Leg {
public:
    /// Throws InputError if the cashflows are not ordered by payment date. An empty leg is allowed
    /// (a fully settled leg) and prices to zero.
    Leg(std::vector<Cashflow> cashflows, Side side, core::Currency currency);

    const std::vector<Cashflow>& cashflows() const noexcept { return cfs_; }
    Side side() const noexcept { return side_; }
    double sign() const noexcept { return side_ == Side::Receive ? 1.0 : -1.0; }
    core::Currency currency() const noexcept { return ccy_; }
    bool empty() const noexcept { return cfs_.empty(); }
    std::size_t size() const noexcept { return cfs_.size(); }
    core::Date lastPayment() const;    ///< throws InputError on an empty leg
    core::Date lastAccrualEnd() const; ///< last accrual end over the coupons (payment for notionals); throws on empty

private:
    std::vector<Cashflow> cfs_;
    Side side_;
    core::Currency ccy_;
};

} // namespace rl::instruments
