#include "rl/models/cheyette.hpp"

#include "rl/core/error.hpp"

#include <algorithm>
#include <cmath>

namespace rl::models {

Cheyette::Cheyette(double kappa, double skew, std::vector<double> times, std::vector<double> levels)
    : kappa_(kappa), skew_(skew), times_(std::move(times)), levels_(std::move(levels)) {
    RL_REQUIRE(std::isfinite(kappa_) && std::isfinite(skew_), "Cheyette: kappa and skew must be finite");
    RL_REQUIRE(!times_.empty() && times_.size() == levels_.size(), "Cheyette: need matching, non-empty level nodes");
    for (std::size_t i = 0; i < times_.size(); ++i) {
        RL_REQUIRE(std::isfinite(times_[i]) && times_[i] > (i ? times_[i - 1] : 0.0), "Cheyette: times must be positive and strictly increasing");
        RL_REQUIRE(std::isfinite(levels_[i]) && levels_[i] > 0.0, "Cheyette: volatility levels must be positive");
    }
}

double Cheyette::level(double t) const noexcept {
    const std::size_t i = static_cast<std::size_t>(std::lower_bound(times_.begin(), times_.end(), t - 1e-12) - times_.begin());
    return levels_[std::min(i, levels_.size() - 1)];
}

double Cheyette::shape(double x) const noexcept { return std::clamp(1.0 + skew_ * x, 0.05, 20.0); }

double Cheyette::G(double tau) const noexcept {
    const double a = kappa_ * tau;
    if (std::fabs(a) < 1e-8) return tau * (1.0 - 0.5 * a);
    return -std::expm1(-a) / kappa_;
}

double Cheyette::bondRatio(double t, double T, double x, double y) const noexcept {
    const double g = G(T - t);
    return std::exp(-g * x - 0.5 * g * g * y);
}

} // namespace rl::models
