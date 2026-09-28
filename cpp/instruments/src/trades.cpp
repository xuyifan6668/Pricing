#include "rl/core/error.hpp"
#include "rl/instruments/basis_swap.hpp"
#include "rl/instruments/bermudan.hpp"
#include "rl/instruments/deposit.hpp"
#include "rl/instruments/fra.hpp"
#include "rl/instruments/future.hpp"
#include "rl/instruments/leg.hpp"
#include "rl/instruments/swap.hpp"
#include "rl/instruments/swaption.hpp"

#include <algorithm>
#include <cmath>

namespace rl::instruments {

namespace {
core::Date firstAccrualStart(const Leg& leg) {
    core::Date first = core::Date::maxDate();
    for (const auto& cf : leg.cashflows())
        std::visit(
            [&](const auto& c) {
                if constexpr (!std::is_same_v<std::decay_t<decltype(c)>, NotionalExchange>) first = std::min(first, c.accrualStart);
            },
            cf);
    return first;
}
core::Date lastAccrualStart(const Leg& leg) {
    core::Date last = core::Date::minDate();
    for (const auto& cf : leg.cashflows())
        std::visit(
            [&](const auto& c) {
                if constexpr (!std::is_same_v<std::decay_t<decltype(c)>, NotionalExchange>) last = std::max(last, c.accrualStart);
            },
            cf);
    return last;
}
} // namespace

// --- Leg -------------------------------------------------------------------------------------------

Leg::Leg(std::vector<Cashflow> cashflows, Side side, core::Currency currency) : cfs_(std::move(cashflows)), side_(side), ccy_(currency) {
    for (std::size_t i = 1; i < cfs_.size(); ++i)
        RL_REQUIRE(paymentDate(cfs_[i]) >= paymentDate(cfs_[i - 1]), "Leg: cashflows must be ordered by payment date");
}

core::Date Leg::lastPayment() const {
    RL_REQUIRE(!cfs_.empty(), "Leg::lastPayment on an empty leg");
    return paymentDate(cfs_.back());
}

core::Date Leg::lastAccrualEnd() const {
    RL_REQUIRE(!cfs_.empty(), "Leg::lastAccrualEnd on an empty leg");
    core::Date last = core::Date::minDate();
    for (const auto& cf : cfs_) last = std::max(last, accrualEndOrPayment(cf));
    return last;
}

// --- Swap ------------------------------------------------------------------------------------------

Swap::Swap(std::string id, Leg fixed, Leg floating) : id_(std::move(id)), fixed_(std::move(fixed)), float_(std::move(floating)) {
    RL_REQUIRE(!id_.empty(), "Swap: id is required");
    RL_REQUIRE(fixed_.side() != float_.side(), "Swap " + id_ + ": legs must have opposite sides");
    RL_REQUIRE(fixed_.currency() == float_.currency(), "Swap " + id_ + ": legs must share one currency (cross-currency is v2)");
    for (const auto& cf : fixed_.cashflows())
        RL_REQUIRE(!std::holds_alternative<OvernightCoupon>(cf) && !std::holds_alternative<TermCoupon>(cf), "Swap " + id_ + ": fixed leg holds a floating coupon");
}

core::Date Swap::maturity() const {
    core::Date m = core::Date::minDate();
    if (!fixed_.empty()) m = std::max(m, fixed_.lastPayment());
    if (!float_.empty()) m = std::max(m, float_.lastPayment());
    return m;
}

// --- BasisSwap -------------------------------------------------------------------------------------

BasisSwap::BasisSwap(std::string id, Leg first, Leg second) : id_(std::move(id)), legs_{std::move(first), std::move(second)} {
    RL_REQUIRE(!id_.empty(), "BasisSwap: id is required");
    RL_REQUIRE(legs_[0].side() != legs_[1].side(), "BasisSwap " + id_ + ": legs must have opposite sides");
    RL_REQUIRE(legs_[0].currency() == legs_[1].currency(), "BasisSwap " + id_ + ": legs must share one currency");
}

core::Date BasisSwap::maturity() const {
    core::Date m = core::Date::minDate();
    for (const auto& l : legs_)
        if (!l.empty()) m = std::max(m, l.lastPayment());
    return m;
}

// --- Deposit / Fra / Future ------------------------------------------------------------------------

Deposit::Deposit(std::string id, core::Index index, core::Date start, core::Date end, double tau, double notional, double rate)
    : id_(std::move(id)), index_(std::move(index)), start_(start), end_(end), tau_(tau), notional_(notional), rate_(rate) {
    RL_REQUIRE(!id_.empty(), "Deposit: id is required");
    RL_REQUIRE(end_ > start_, "Deposit " + id_ + ": end must follow start");
    RL_REQUIRE(tau_ > 0.0 && std::isfinite(rate_), "Deposit " + id_ + ": bad tau or rate");
}

Fra::Fra(std::string id, core::Index index, core::Date fixingDate, core::Date start, core::Date end, double tau, double notional, double rate)
    : id_(std::move(id)), index_(std::move(index)), fixing_(fixingDate), start_(start), end_(end), tau_(tau), notional_(notional), rate_(rate) {
    RL_REQUIRE(!id_.empty(), "Fra: id is required");
    RL_REQUIRE(end_ > start_ && fixing_ <= start_, "Fra " + id_ + ": need fixing <= start < end");
    RL_REQUIRE(tau_ > 0.0 && std::isfinite(rate_), "Fra " + id_ + ": bad tau or rate");
}

Future::Future(std::string id, core::Index index, core::Date start, core::Date end, double tau, double price, double convexityAdjustment)
    : id_(std::move(id)), index_(std::move(index)), start_(start), end_(end), tau_(tau), price_(price), convexity_(convexityAdjustment) {
    RL_REQUIRE(!id_.empty(), "Future: id is required");
    RL_REQUIRE(end_ > start_ && tau_ > 0.0, "Future " + id_ + ": bad reference period");
    RL_REQUIRE(price_ > 0.0 && price_ < 200.0 && std::isfinite(convexity_), "Future " + id_ + ": bad price or convexity");
}

// --- Swaption / Bermudan ---------------------------------------------------------------------------

Swaption::Swaption(std::string id, Swap underlying, core::Date expiry, Settlement settlement)
    : id_(std::move(id)), swap_(std::move(underlying)), expiry_(expiry), settlement_(settlement) {
    RL_REQUIRE(!id_.empty(), "Swaption: id is required");
    RL_REQUIRE(expiry_ <= firstAccrualStart(swap_.fixedLeg()), "Swaption " + id_ + ": expiry after the underlying starts");
}

BermudanSwaption::BermudanSwaption(std::string id, Swap underlying, std::vector<core::Date> exerciseDates, int noticeDays)
    : id_(std::move(id)), swap_(std::move(underlying)), exercise_(std::move(exerciseDates)), notice_(noticeDays) {
    RL_REQUIRE(!id_.empty(), "BermudanSwaption: id is required");
    RL_REQUIRE(!exercise_.empty(), "BermudanSwaption " + id_ + ": no exercise dates");
    RL_REQUIRE(notice_ >= 0, "BermudanSwaption " + id_ + ": notice must be non-negative");
    for (std::size_t i = 1; i < exercise_.size(); ++i) RL_REQUIRE(exercise_[i] > exercise_[i - 1], "BermudanSwaption " + id_ + ": exercise dates not increasing");
    RL_REQUIRE(exercise_.back() <= lastAccrualStart(swap_.fixedLeg()), "BermudanSwaption " + id_ + ": exercise after the last period start");
}

} // namespace rl::instruments
