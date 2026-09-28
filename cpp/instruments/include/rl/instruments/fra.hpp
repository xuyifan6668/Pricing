#pragma once
#include "rl/core/date.hpp"
#include "rl/core/index.hpp"

#include <string>

namespace rl::instruments {

/// Forward rate agreement on a term index (or a single-period forward on an overnight index):
/// pays notional · tau · (fixing − rate) at the end of the accrual period.
/// SIMPLIFICATION: settled at accrual end, not the discounted-at-start FRA convention; irrelevant for
/// curve construction, where only the implied forward matters.
class Fra {
public:
    Fra(std::string id, core::Index index, core::Date fixingDate, core::Date start, core::Date end, double tau, double notional, double rate);

    const std::string& id() const noexcept { return id_; }
    const core::Index& index() const noexcept { return index_; }
    core::Date fixingDate() const noexcept { return fixing_; }
    core::Date start() const noexcept { return start_; }
    core::Date end() const noexcept { return end_; }
    double tau() const noexcept { return tau_; }
    double notional() const noexcept { return notional_; }
    double rate() const noexcept { return rate_; }

private:
    std::string id_;
    core::Index index_;
    core::Date fixing_, start_, end_;
    double tau_, notional_, rate_;
};

} // namespace rl::instruments
