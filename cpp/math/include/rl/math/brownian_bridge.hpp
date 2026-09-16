#pragma once
// Brownian bridge construction of a discretely observed Brownian path from i.i.d. N(0,1) draws,
// in the order that puts the most variance into the first coordinates (terminal value first, then
// midpoints, ...). With Sobol, the low dimensions are the well-distributed ones, so this is what
// turns a Sobol sequence into a genuinely better path than the incremental construction.
//
// Reference: Glasserman (2004) §3.1; Jäckel (2002) ch. 10 (the index/weight precomputation below is
// his). Everything here is passive: the draws are doubles, the path is doubles; the model's
// x_t = f(W_t; θ) is where T enters.

#include <cstddef>
#include <span>
#include <vector>

namespace rl::math {

class BrownianBridge {
public:
    /// `times` strictly increasing, all > 0 (t₀ = 0 with W₀ = 0 is implicit). Size = number of steps.
    explicit BrownianBridge(std::vector<double> times);

    std::size_t size() const noexcept { return t_.size(); }
    const std::vector<double>& times() const noexcept { return t_; }

    /// gaussians (size n, i.i.d. N(0,1)) → path W(t₀)..W(tₙ₋₁) (size n). Coordinate 0 of `gaussians`
    /// drives W(tₙ₋₁), coordinate 1 the midpoint, and so on.
    void buildPath(std::span<const double> gaussians, std::span<double> path) const;
    /// Same, but returns the increments ΔWᵢ = W(tᵢ) − W(tᵢ₋₁) (what an Euler/exact-step model consumes).
    void buildIncrements(std::span<const double> gaussians, std::span<double> increments) const;

    // Precomputed tables (exposed for tests): for bridge step i (i = 1..n−1), fill point index_[i]
    // between leftIndex_[i] and rightIndex_[i] with weights and stdDev.
    const std::vector<std::size_t>& bridgeIndex() const noexcept { return index_; }
    const std::vector<std::size_t>& leftIndex() const noexcept { return leftIndex_; }
    const std::vector<std::size_t>& rightIndex() const noexcept { return rightIndex_; }
    const std::vector<double>& leftWeight() const noexcept { return leftWeight_; }
    const std::vector<double>& rightWeight() const noexcept { return rightWeight_; }
    const std::vector<double>& stdDev() const noexcept { return stdDev_; }

private:
    void build(); // the maths
    std::vector<double> t_;
    std::vector<std::size_t> index_, leftIndex_, rightIndex_;
    std::vector<double> leftWeight_, rightWeight_, stdDev_;
};

} // namespace rl::math
