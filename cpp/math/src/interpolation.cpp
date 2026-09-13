#include "rl/math/detail/interpolation.ipp"

#include "rl/core/error.hpp"

#include <algorithm>

namespace rl::math::detail {

std::size_t segment(const std::vector<double>& xs, double x) noexcept {
    // upper_bound gives the first node > x; the segment starts one before it.
    auto it = std::upper_bound(xs.begin(), xs.end(), x);
    if (it == xs.begin()) return 0;
    std::size_t i = static_cast<std::size_t>(it - xs.begin()) - 1;
    return std::min(i, xs.size() - 2);
}

void validateNodes(const std::vector<double>& xs, std::size_t ny, std::size_t minPoints) {
    RL_REQUIRE(xs.size() >= minPoints, "interpolation needs at least " + std::to_string(minPoints) + " nodes");
    RL_REQUIRE(xs.size() == ny, "interpolation: xs and ys differ in size");
    for (std::size_t i = 1; i < xs.size(); ++i)
        RL_REQUIRE(xs[i] > xs[i - 1], "interpolation nodes must be strictly increasing");
    for (double x : xs) RL_REQUIRE(std::isfinite(x), "interpolation node is not finite");
}

} // namespace rl::math::detail

RL_MATH_INSTANTIATE_INTERPOLATION(double)
