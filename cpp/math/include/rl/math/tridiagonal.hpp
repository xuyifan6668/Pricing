#pragma once
// Tridiagonal linear systems. The Crank–Nicolson step  is one solve per time step, and
// tapes through it (T = Real), so the solve must be a plain loop over T — no Eigen here.

#include "rl/math/scalar.hpp"

#include <cstddef>
#include <vector>

namespace rl::math {

/// A x = d with A tridiagonal: sub-diagonal a (size n−1), diagonal b (size n), super-diagonal c (size n−1).
///   b₀ x₀ + c₀ x₁                  = d₀
///   aᵢ₋₁ xᵢ₋₁ + bᵢ xᵢ + cᵢ xᵢ₊₁   = dᵢ      0 < i < n−1
///   aₙ₋₂ xₙ₋₂ + bₙ₋₁ xₙ₋₁          = dₙ₋₁
template <Scalar T>
class Tridiagonal {
public:
    Tridiagonal(std::vector<T> sub, std::vector<T> diag, std::vector<T> super);

    std::size_t size() const noexcept { return diag_.size(); }
    const std::vector<T>& sub() const noexcept { return sub_; }
    const std::vector<T>& diag() const noexcept { return diag_; }
    const std::vector<T>& super() const noexcept { return super_; }

    /// True if |bᵢ| ≥ |aᵢ₋₁| + |cᵢ| for every row (strict for at least one) — the Thomas algorithm is
    /// then stable without pivoting. The PDE operator satisfies this for every sane Δt/Δx.
    bool isDiagonallyDominant() const;

    /// Solve A x = d. Throws InputError on size mismatch, NumericalError on a zero pivot.
    std::vector<T> solve(const std::vector<T>& d) const;
    /// In-place variant for the hot loop: x is overwritten; `work` must have size() entries.
    void solveInPlace(std::vector<T>& x, std::vector<T>& work) const;

    /// y = A x
    std::vector<T> multiply(const std::vector<T>& x) const;

private:
    std::vector<T> sub_, diag_, super_;
};

} // namespace rl::math
