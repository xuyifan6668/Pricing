#pragma once
// Gauss–Hermite quadrature, normalised for expectations under N(0,1):
//   E[g(Z)] ≈ Σᵢ wᵢ g(xᵢ)
// Nodes/weights are passive doubles (they depend only on n). The integrand is a template so the
// expectation is a T when g returns T: that is how the LGM swaption cross-check and the martingale
// test  run under both scalars.

#include "rl/math/scalar.hpp"

#include <cstddef>
#include <vector>

namespace rl::math {

class GaussHermite {
public:
    /// n-point rule (n ≥ 1). Exact for polynomials of degree ≤ 2n − 1 against the Gaussian weight.
    explicit GaussHermite(std::size_t n);

    std::size_t size() const noexcept { return x_.size(); }
    const std::vector<double>& nodes() const noexcept { return x_; }     ///< xᵢ, ascending, symmetric about 0
    const std::vector<double>& weights() const noexcept { return w_; }   ///< wᵢ > 0, Σ wᵢ = 1

    /// Σ wᵢ g(xᵢ); g: double → T.
    template <class G>
    auto integrate(G&& g) const -> decltype(g(0.0)) {
        using R = decltype(g(0.0));
        R acc = R(0.0);
        for (std::size_t i = 0; i < x_.size(); ++i) acc += w_[i] * g(x_[i]);
        return acc;
    }

    /// E[g(μ + σ Z)] — the form the LGM code needs (Z ~ N(0,1), x_t = μ + σ Z).
    template <class G>
    auto integrate(double mu, double sigma, G&& g) const -> decltype(g(0.0)) {
        return integrate([&](double z) { return g(mu + sigma * z); });
    }

private:
    void build(std::size_t n); // the maths; see quadrature.cpp
    std::vector<double> x_, w_;
};

} // namespace rl::math
