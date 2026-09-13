#pragma once
#include "rl/core/error.hpp"
#include "rl/math/tridiagonal.hpp"

#include <cmath>
#include <string>
#include <utility>

namespace rl::math {

template <Scalar T>
Tridiagonal<T>::Tridiagonal(std::vector<T> sub, std::vector<T> diag, std::vector<T> super)
    : sub_(std::move(sub)), diag_(std::move(diag)), super_(std::move(super)) {
    RL_REQUIRE(diag_.size() >= 1, "tridiagonal system must have at least one row");
    RL_REQUIRE(sub_.size() == diag_.size() - 1 && super_.size() == diag_.size() - 1, "tridiagonal band sizes must be n-1, n, n-1");
}

/**
 * Row-wise check |bᵢ| ≥ |aᵢ₋₁| + |cᵢ| (weak), with strict inequality somewhere. Use value().
 */
template <Scalar T>
bool Tridiagonal<T>::isDiagonallyDominant() const {
    const std::size_t n = size();
    bool strict = false;
    for (std::size_t i = 0; i < n; ++i) {
        double off = 0.0;
        if (i > 0) off += std::fabs(value(sub_[i - 1]));
        if (i + 1 < n) off += std::fabs(value(super_[i]));
        const double d = std::fabs(value(diag_[i]));
        if (d < off) return false;
        if (d > off) strict = true;
    }
    return strict;
}

/**
 * Thomas algorithm (Gaussian elimination without pivoting specialised to the band). Forward sweep:
 *   c'₀ = c₀/b₀,  d'₀ = d₀/b₀
 *   c'ᵢ = cᵢ / (bᵢ − aᵢ₋₁ c'ᵢ₋₁),   d'ᵢ = (dᵢ − aᵢ₋₁ d'ᵢ₋₁) / (bᵢ − aᵢ₋₁ c'ᵢ₋₁)
 * Back substitution:
 *   xₙ₋₁ = d'ₙ₋₁,   xᵢ = d'ᵢ − c'ᵢ xᵢ₊₁
 * 8n flops, O(n) memory for the modified coefficients (that is `work`). Throw NumericalError if a
 * pivot bᵢ − aᵢ₋₁ c'ᵢ₋₁ has |value| < 1e-300 (only possible if the matrix is singular or far from
 * dominant). This is THE hot loop of the PDE pricer (one call per time step, n ≈ 500), and the AAD scalar
 * tapes it: no branching on T, no allocation inside.
 * Test: against Eigen's dense solve on random dominant systems to 1e-12; multiply(solve(d)) == d.
 */
template <Scalar T>
void Tridiagonal<T>::solveInPlace(std::vector<T>& x, std::vector<T>& work) const {
    // Thomas algorithm, Numerical Recipes 3e §2.4 (tridag). `work` holds the
    // modified super-diagonal. Pivots are checked on value() so a Real solve records no branch.
    const std::size_t n = size();
    RL_REQUIRE(x.size() == n && work.size() == n, "tridiagonal solveInPlace: size mismatch");
    T piv = diag_[0];
    if (std::fabs(value(piv)) < 1e-300) throw NumericalError("tridiagonal solve: zero pivot at row 0");
    x[0] = x[0] / piv;
    for (std::size_t i = 1; i < n; ++i) {
        work[i] = super_[i - 1] / piv;
        piv = diag_[i] - sub_[i - 1] * work[i];
        if (std::fabs(value(piv)) < 1e-300) throw NumericalError("tridiagonal solve: zero pivot at row " + std::to_string(i));
        x[i] = (x[i] - sub_[i - 1] * x[i - 1]) / piv;
    }
    for (std::size_t i = n - 1; i-- > 0;) x[i] -= work[i + 1] * x[i + 1];
}

template <Scalar T>
std::vector<T> Tridiagonal<T>::solve(const std::vector<T>& d) const {
    RL_REQUIRE(d.size() == size(), "tridiagonal solve: rhs size mismatch");
    std::vector<T> x(d);
    std::vector<T> work(size());
    solveInPlace(x, work);
    return x;
}

/**
 * yᵢ = aᵢ₋₁ xᵢ₋₁ + bᵢ xᵢ + cᵢ xᵢ₊₁ (with the missing terms at the ends dropped). Plain plumbing but
 * it is the explicit half of Crank–Nicolson, so keep it allocation-free apart from the result.
 */
template <Scalar T>
std::vector<T> Tridiagonal<T>::multiply(const std::vector<T>& x) const {
    RL_REQUIRE(x.size() == size(), "tridiagonal multiply: size mismatch");
    const std::size_t n = size();
    std::vector<T> y(n);
    for (std::size_t i = 0; i < n; ++i) {
        T acc = diag_[i] * x[i];
        if (i > 0) acc += sub_[i - 1] * x[i - 1];
        if (i + 1 < n) acc += super_[i] * x[i + 1];
        y[i] = acc;
    }
    return y;
}

} // namespace rl::math

#define RL_MATH_INSTANTIATE_TRIDIAGONAL(T) template class rl::math::Tridiagonal<T>;
