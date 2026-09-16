#include "rl/math/brownian_bridge.hpp"

#include "rl/core/error.hpp"

#include <cmath>
#include <utility>

namespace rl::math {

BrownianBridge::BrownianBridge(std::vector<double> times) : t_(std::move(times)) {
    RL_REQUIRE(!t_.empty(), "Brownian bridge needs at least one time");
    RL_REQUIRE(t_.front() > 0.0, "Brownian bridge times must be > 0");
    for (std::size_t i = 1; i < t_.size(); ++i) RL_REQUIRE(t_[i] > t_[i - 1], "Brownian bridge times must be strictly increasing");
    const std::size_t n = t_.size();
    index_.assign(n, 0);
    leftIndex_.assign(n, 0);
    rightIndex_.assign(n, 0);
    leftWeight_.assign(n, 0.0);
    rightWeight_.assign(n, 0.0);
    stdDev_.assign(n, 0.0);
    build();
}

// Jäckel (2002) §10.8.3 / Glasserman (2004) §3.1 bridge tables.
void BrownianBridge::build() {
    const std::size_t n = t_.size();
    std::vector<std::size_t> map(n, 0); // map[k] = 1 + bridge step that fills slot k; 0 = unfilled
    const std::size_t origin = n;       // sentinel for W(0) = 0
    map[n - 1] = 1;
    index_[0] = n - 1;
    leftIndex_[0] = origin;
    rightIndex_[0] = n - 1;
    leftWeight_[0] = 0.0;
    rightWeight_[0] = 0.0;
    stdDev_[0] = std::sqrt(t_[n - 1]);
    std::size_t j = 0;
    for (std::size_t i = 1; i < n; ++i) {
        while (map[j] != 0) ++j;          // next unfilled slot
        std::size_t k = j;
        while (map[k] == 0) ++k;          // next filled slot after it
        const std::size_t l = j + (k - 1 - j) / 2; // midpoint slot
        map[l] = i + 1;
        const std::size_t left = (j == 0) ? origin : j - 1;
        const double ta = (j == 0) ? 0.0 : t_[j - 1], tl = t_[l], tk = t_[k];
        index_[i] = l;
        leftIndex_[i] = left;
        rightIndex_[i] = k;
        leftWeight_[i] = (tk - tl) / (tk - ta);
        rightWeight_[i] = (tl - ta) / (tk - ta);
        stdDev_[i] = std::sqrt((tl - ta) * (tk - tl) / (tk - ta));
        j = k + 1;
        if (j >= n) j = 0;
    }
}

void BrownianBridge::buildPath(std::span<const double> gaussians, std::span<double> path) const {
    const std::size_t n = t_.size();
    RL_REQUIRE(gaussians.size() >= n && path.size() >= n, "BrownianBridge::buildPath: span too small");
    path[index_[0]] = stdDev_[0] * gaussians[0];
    for (std::size_t i = 1; i < n; ++i) {
        const double left = leftIndex_[i] == n ? 0.0 : path[leftIndex_[i]];
        path[index_[i]] = leftWeight_[i] * left + rightWeight_[i] * path[rightIndex_[i]] + stdDev_[i] * gaussians[i];
    }
}

void BrownianBridge::buildIncrements(std::span<const double> gaussians, std::span<double> increments) const {
    buildPath(gaussians, increments);
    for (std::size_t i = t_.size(); i-- > 1;) increments[i] -= increments[i - 1];
}

} // namespace rl::math
