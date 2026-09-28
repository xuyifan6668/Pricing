#pragma once
// Bermudan swaption: id + underlying swap (owned) + exercise dates + notice. Exercise on date e_i
// enters the underlying's coupons whose accrual starts on or after the exercise's effective date;
// the LgmBermudanProduct<T> does the slicing (remainingCouponsFrom() below is the helper it uses).
//
// SIMPLIFICATION: v1 treats notice as a whole number of business days on the swap's calendar and
// exercises into the swap "as is" (no fee, no partial exercise), which is the vanilla USD Bermudan.

#include "rl/core/date.hpp"
#include "rl/instruments/swap.hpp"

#include <string>
#include <vector>

namespace rl::instruments {

class BermudanSwaption {
public:
    /// Throws InputError unless exercise dates are strictly increasing, non-empty, and each is before
    /// the underlying's last accrual start (an exercise after the last period start is meaningless).
    BermudanSwaption(std::string id, Swap underlying, std::vector<core::Date> exerciseDates, int noticeDays);

    const std::string& id() const noexcept { return id_; }
    const Swap& underlying() const noexcept { return swap_; }
    const std::vector<core::Date>& exerciseDates() const noexcept { return exercise_; }
    int noticeDays() const noexcept { return notice_; }
    bool payer() const noexcept { return swap_.payerOfFixed(); }

private:
    std::string id_;
    Swap swap_;
    std::vector<core::Date> exercise_;
    int notice_;
};

} // namespace rl::instruments
