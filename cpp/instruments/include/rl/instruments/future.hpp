#pragma once
#include "rl/core/date.hpp"
#include "rl/core/index.hpp"

#include <string>

namespace rl::instruments {

/// Three-month compounded-overnight future (CME SR3): price = 100·(1 − R), R the daily-compounded rate
/// over the reference period [start, end] (IMM date to next IMM date), ACT/360.
/// The convexity adjustment c is a *contract term here* (rate units): futures rate = forward + c.
/// SIMPLIFICATION: production derives c from a short-rate model (Hull–White: grows with expiry and
/// the reference period); v1 takes it from the quote file.
class Future {
public:
    Future(std::string id, core::Index index, core::Date start, core::Date end, double tau, double price, double convexityAdjustment);

    const std::string& id() const noexcept { return id_; }
    const core::Index& index() const noexcept { return index_; }
    core::Date start() const noexcept { return start_; }
    core::Date end() const noexcept { return end_; }
    double tau() const noexcept { return tau_; }
    double price() const noexcept { return price_; }
    double convexityAdjustment() const noexcept { return convexity_; }
    /// 1 − price/100: the quoted futures rate.
    double quotedRate() const noexcept { return 1.0 - price_ / 100.0; }

private:
    std::string id_;
    core::Index index_;
    core::Date start_, end_;
    double tau_, price_, convexity_;
};

} // namespace rl::instruments
