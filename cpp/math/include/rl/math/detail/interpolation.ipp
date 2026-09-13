#pragma once
// Template definitions for interpolation.hpp. Included ONLY by the explicit-instantiation TUs
// (cpp/math/src/interpolation.cpp for double, cpp/aad/src/instantiate_math.cpp for Real).

#include "rl/core/error.hpp"
#include "rl/math/interpolation.hpp"
#include "rl/math/tridiagonal.hpp"

#include <algorithm>
#include <cmath>
#include <utility>

namespace rl::math {


// ------------------------------------------------------------------------------------------------
// Linear
// ------------------------------------------------------------------------------------------------

template <Scalar T>
Linear<T>::Linear(std::vector<double> xs, std::vector<T> ys, Extrapolation extrap)
    : xs_(std::move(xs)), ys_(std::move(ys)), extrap_(extrap) {
    detail::validateNodes(xs_, ys_.size(), 2);
}

/**
 * Piecewise-linear value. With i = segment(x):
 *   y(x) = y_i + (y_{i+1} − y_i) · (x − x_i)/(x_{i+1} − x_i)
 * Extrapolation::Flat clamps x to [x_0, x_{n−1}] first; Extrapolation::Linear uses the end segment
 * as is (the formula above already does that since segment() clamps i).
 * Hot in the bootstrap (called per residual per node); O(log n) from the bisection, no allocation.
 * Test: reproduces nodes exactly; midpoint is the average; both extrapolation modes.
 */
template <Scalar T>
T Linear<T>::operator()(double x) const {
    if (extrap_ == Extrapolation::Flat) {
        if (x <= xs_.front()) return ys_.front();
        if (x >= xs_.back()) return ys_.back();
    }
    const std::size_t i = detail::segment(xs_, x);
    const double w = (x - xs_[i]) / (xs_[i + 1] - xs_[i]);
    return ys_[i] + (ys_[i + 1] - ys_[i]) * w;
}

/**
 * dy/dx on the segment containing x: (y_{i+1} − y_i)/(x_{i+1} − x_i). Zero outside the range under
 * Flat extrapolation. At a node, use the right-hand segment (the convention downstream relies on:
 * forward(t_i) is the forward for (t_i, t_{i+1}]).
 */
template <Scalar T>
T Linear<T>::derivative(double x) const {
    if (extrap_ == Extrapolation::Flat && (x <= xs_.front() || x >= xs_.back())) return T(0.0);
    const std::size_t i = detail::segment(xs_, x);
    return (ys_[i + 1] - ys_[i]) / (xs_[i + 1] - xs_[i]);
}

// ------------------------------------------------------------------------------------------------
// LogLinearDf
// ------------------------------------------------------------------------------------------------

namespace detail {
template <Scalar T>
std::vector<T> logAll(const std::vector<T>& dfs) {
    std::vector<T> out;
    out.reserve(dfs.size());
    for (const auto& p : dfs) {
        RL_REQUIRE(value(p) > 0.0, "discount factors must be positive for log-linear interpolation");
        out.push_back(log(p));
    }
    return out;
}
} // namespace detail

template <Scalar T>
LogLinearDf<T>::LogLinearDf(std::vector<double> xs, std::vector<T> dfs, Extrapolation extrap)
    : logInterp_(std::move(xs), detail::logAll(dfs), extrap) {}

/** P(x) = exp(linear interpolation of log P): nodes reproduced, constant forward between nodes. */
template <Scalar T>
T LogLinearDf<T>::operator()(double x) const {
    return exp(logInterp_(x));
}

/** dP/dx = P(x) · d(log P)/dx = −f(x) P(x). */
template <Scalar T>
T LogLinearDf<T>::derivative(double x) const {
    return logInterp_.derivative(x) * exp(logInterp_(x));
}

/** f(x) = −d(log P)/dx, piecewise constant. */
template <Scalar T>
T LogLinearDf<T>::forward(double x) const {
    return -logInterp_.derivative(x);
}

// ------------------------------------------------------------------------------------------------
// MonotoneConvex — Hagan & West (2006)
// ------------------------------------------------------------------------------------------------

template <Scalar T>
MonotoneConvex<T>::MonotoneConvex(std::vector<double> ts, std::vector<T> ys, bool enforcePositivity)
    : ts_(std::move(ts)), ys_(std::move(ys)), positivity_(enforcePositivity) {
    detail::validateNodes(ts_, ys_.size(), 2);
    RL_REQUIRE(ts_.front() == 0.0 && value(ys_.front()) == 0.0, "MonotoneConvex needs the first node at (t=0, y=0)");
    build();
}

namespace detail {

// The Hagan–West correction g(x) on one interval and its integral G(x) = ∫_0^x g, x ∈ [0, 1].
// g0 = f_{i−1} − f^d_i and g1 = f_i − f^d_i are the deviations of the two node forwards from the
// interval's discrete forward. Every regime satisfies g(0) = g0, g(1) = g1 and G(1) = 0; the last
// property is what makes the interpolant reproduce the node ordinates exactly.
//
// The regime is chosen on value(g0), value(g1): it is an active-set decision, so it is frozen on
// the tape exactly like a max(). The four regimes (Hagan & West 2006, §5.2):
//   (i)   g1 between −½g0 and −2g0: a plain quadratic is already monotone.
//   (ii)  g1 beyond −2g0: flat at g0 up to η, then a quadratic up to g1.
//   (iii) g1 between 0 and −½g0: quadratic down to g1 at η, then flat.
//   (iv)  g0, g1 of the same sign: two quadratics meeting at a common level A at η.
template <Scalar T>
struct MonotoneConvexPiece {
    T g;
    T G;
};

template <Scalar T>
MonotoneConvexPiece<T> monotoneConvexPiece(const T& g0, const T& g1, double x) {
    // Deviations below 1e-13 (in rate units) are rounding noise from two equal discrete forwards.
    // DESIGN: with one deviation exactly zero, regimes (iii)/(iv) degenerate to η → 0 or 1, i.e. a
    // jump in the forward at the node. The plain quadratic is used there instead: it still meets
    // both end values and integrates to zero, and it keeps the forward continuous where two
    // neighbouring intervals happen to have the same average forward (a flat long end).
    constexpr double kNoise = 1e-13;
    const double a = std::abs(value(g0)) < kNoise ? 0.0 : value(g0);
    const double b = std::abs(value(g1)) < kNoise ? 0.0 : value(g1);
    const double x2 = x * x, x3 = x2 * x;

    // Both zero (a locally flat curve) also takes the quadratic: its value is zero there, and being
    // linear in (g0, g1) it keeps the tape's dependence on the nodes, so a Newton step from a flat
    // initial guess still sees how the forward responds to each node.
    const bool quadratic = a == 0.0 || b == 0.0 || (a < 0.0 && -0.5 * a <= b && b <= -2.0 * a) ||
                           (a > 0.0 && -0.5 * a >= b && b >= -2.0 * a);
    if (quadratic) {
        return {g0 * (1.0 - 4.0 * x + 3.0 * x2) + g1 * (-2.0 * x + 3.0 * x2),
                g0 * (x - 2.0 * x2 + x3) + g1 * (-x2 + x3)};
    }
    // In regimes (ii)–(iv) the breakpoint η is itself a function of (g0, g1), so it is carried as a T:
    // dropping its sensitivity would make the curve's node derivatives wrong while leaving every
    // value right. Only the *comparison* of x against η is done on values.
    // (ii)
    if ((a < 0.0 && b > -2.0 * a) || (a > 0.0 && b < -2.0 * a)) {
        const T eta = (g1 + 2.0 * g0) / (g1 - g0);
        if (x <= value(eta)) return {g0, g0 * x};
        const T u = (x - eta) / (1.0 - eta);
        return {g0 + (g1 - g0) * (u * u), g0 * x + (g1 - g0) * ((1.0 - eta) / 3.0 * (u * u * u))};
    }
    // (iii)
    if ((a > 0.0 && b < 0.0 && b > -0.5 * a) || (a < 0.0 && b > 0.0 && b < -0.5 * a)) {
        const T eta = 3.0 * g1 / (g1 - g0);
        if (x >= value(eta)) return {g1, g1 * x + (g0 - g1) * (eta / 3.0)};
        const T u = (eta - x) / eta;
        return {g1 + (g0 - g1) * (u * u), g1 * x + (g0 - g1) * (eta / 3.0 * (1.0 - u * u * u))};
    }
    // (iv): same sign.
    const T eta = g1 / (g0 + g1);
    const T A = -1.0 * g0 * g1 / (g0 + g1);
    if (x <= value(eta)) {
        const T u = (eta - x) / eta;
        return {A + (g0 - A) * (u * u), A * x + (g0 - A) * (eta / 3.0 * (1.0 - u * u * u))};
    }
    const T u = (x - eta) / (1.0 - eta);
    return {A + (g1 - A) * (u * u), A * x + (g0 - A) * (eta / 3.0) + (g1 - A) * ((1.0 - eta) / 3.0 * (u * u * u))};
}

/// min(max(x, 0), hi) by value comparison, keeping the chosen operand's tape node.
template <Scalar T>
T collar(const T& x, const T& hi) {
    if (value(x) < 0.0) return T(0.0);
    return value(x) > value(hi) ? hi : x;
}

} // namespace detail

/**
 * Discrete forwards and node forwards (Hagan–West §5, eq. (22)–(26)).
 *   f^d_i = (y_i − y_{i−1}) / (t_i − t_{i−1})                                       i = 1..n−1
 *   f_i   = [ (t_i − t_{i−1}) f^d_{i+1} + (t_{i+1} − t_i) f^d_i ] / (t_{i+1} − t_{i−1})   interior
 *   f_0   = f^d_1 − ½ (f_1 − f^d_1),   f_{n−1} = f^d_{n−1} − ½ (f_{n−2} − f^d_{n−1})
 * Positivity collar (§6): f_i ← min(max(f_i, 0), 2 min(f^d_i, f^d_{i+1})); the ends use their one
 * neighbour. The collar is what stops the interpolant dipping below zero between nodes when the
 * curve is steep; on negative-rate curves it must be switched off (enforcePositivity = false).
 */
template <Scalar T>
void MonotoneConvex<T>::build() {
    const std::size_t n = ts_.size();
    fd_.resize(n - 1);
    for (std::size_t i = 1; i < n; ++i) fd_[i - 1] = (ys_[i] - ys_[i - 1]) / (ts_[i] - ts_[i - 1]);

    f_.assign(n, T(0.0));
    for (std::size_t i = 1; i + 1 < n; ++i) {
        const double hl = ts_[i] - ts_[i - 1], hr = ts_[i + 1] - ts_[i];
        f_[i] = (hl * fd_[i] + hr * fd_[i - 1]) / (hl + hr); // fd_[i] is f^d_{i+1}
    }
    if (positivity_) {
        for (std::size_t i = 1; i + 1 < n; ++i)
            f_[i] = detail::collar(f_[i], 2.0 * min(fd_[i - 1], fd_[i]));
    }
    if (n == 2) {
        f_[0] = f_[1] = fd_[0];
        return;
    }
    f_[0] = fd_[0] - 0.5 * (f_[1] - fd_[0]);
    f_[n - 1] = fd_[n - 2] - 0.5 * (f_[n - 2] - fd_[n - 2]);
    if (positivity_) {
        f_[0] = detail::collar(f_[0], 2.0 * fd_[0]);
        f_[n - 1] = detail::collar(f_[n - 1], 2.0 * fd_[n - 2]);
    }
}

/**
 * f(t) = f^d_i + g(x) on (t_{i−1}, t_i], x = (t − t_{i−1})/(t_i − t_{i−1}); flat at f_{n−1} beyond the
 * last node (flat-forward extrapolation). Continuous at the nodes because g(1) on one interval and
 * g(0) on the next both bring the forward to the shared node forward f_i.
 */
template <Scalar T>
T MonotoneConvex<T>::forward(double t) const {
    if (t <= ts_.front()) return f_.front();
    if (t >= ts_.back()) return f_.back();
    const std::size_t i = detail::segment(ts_, t);
    const double x = (t - ts_[i]) / (ts_[i + 1] - ts_[i]);
    const T& fd = fd_[i];
    return fd + detail::monotoneConvexPiece<T>(f_[i] - fd, f_[i + 1] - fd, x).g;
}

/**
 * y(t) = y_{i−1} + (t_i − t_{i−1}) [ f^d_i x + G(x) ], G = ∫_0^x g. G(1) = 0 in every regime, so
 * y(t_i) = y_i exactly and the bootstrap reprices to machine precision under this scheme too.
 */
template <Scalar T>
T MonotoneConvex<T>::operator()(double t) const {
    if (t <= ts_.front()) return ys_.front() + f_.front() * (t - ts_.front());
    if (t >= ts_.back()) return ys_.back() + f_.back() * (t - ts_.back());
    const std::size_t i = detail::segment(ts_, t);
    const double h = ts_[i + 1] - ts_[i];
    const double x = (t - ts_[i]) / h;
    const T& fd = fd_[i];
    return ys_[i] + h * (fd * x + detail::monotoneConvexPiece<T>(f_[i] - fd, f_[i + 1] - fd, x).G);
}

// ------------------------------------------------------------------------------------------------
// NaturalCubic
// ------------------------------------------------------------------------------------------------

template <Scalar T>
NaturalCubic<T>::NaturalCubic(std::vector<double> xs, std::vector<T> ys, Extrapolation extrap)
    : xs_(std::move(xs)), ys_(std::move(ys)), extrap_(extrap) {
    detail::validateNodes(xs_, ys_.size(), 3);
    build();
}

// Numerical Recipes 3e §3.3 / Burden & Faires §3.5. Node second derivatives M_i
// from the (n−2)×(n−2) symmetric diagonally dominant tridiagonal system
//   h_{i−1} M_{i−1} + 2(h_{i−1} + h_i) M_i + h_i M_{i+1} = 6[(y_{i+1} − y_i)/h_i − (y_i − y_{i−1})/h_{i−1}],
// with M_0 = M_{n−1} = 0. Used only by the vol cube; curves use Linear / MonotoneConvex.
template <Scalar T>
void NaturalCubic<T>::build() {
    const std::size_t n = xs_.size();
    m_.assign(n, T(0.0));
    if (n < 3) return;
    const std::size_t k = n - 2;
    std::vector<T> sub(k > 1 ? k - 1 : 0), diag(k), super(k > 1 ? k - 1 : 0), rhs(k);
    for (std::size_t i = 1; i <= k; ++i) {
        const double h0 = xs_[i] - xs_[i - 1], h1 = xs_[i + 1] - xs_[i];
        diag[i - 1] = T(2.0 * (h0 + h1));
        if (i > 1) sub[i - 2] = T(h0);
        if (i < k) super[i - 1] = T(h1);
        rhs[i - 1] = 6.0 * ((ys_[i + 1] - ys_[i]) / h1 - (ys_[i] - ys_[i - 1]) / h0);
    }
    Tridiagonal<T> A(std::move(sub), std::move(diag), std::move(super));
    std::vector<T> m = A.solve(rhs);
    for (std::size_t i = 0; i < k; ++i) m_[i + 1] = m[i];
}

/**
 * On segment i with h = x_{i+1} − x_i, a = (x_{i+1} − x)/h, b = (x − x_i)/h:
 *   y(x) = a y_i + b y_{i+1} + [(a³ − a) M_i + (b³ − b) M_{i+1}] h²/6
 * Flat extrapolation clamps x; Linear continues with derivative(x_end).
 */
template <Scalar T>
T NaturalCubic<T>::operator()(double x) const {
    if (extrap_ == Extrapolation::Flat) {
        if (x <= xs_.front()) return ys_.front();
        if (x >= xs_.back()) return ys_.back();
    } else if (x < xs_.front()) {
        return ys_.front() + derivative(xs_.front()) * (x - xs_.front());
    } else if (x > xs_.back()) {
        return ys_.back() + derivative(xs_.back()) * (x - xs_.back());
    }
    const std::size_t i = detail::segment(xs_, x);
    const double h = xs_[i + 1] - xs_[i];
    const double a = (xs_[i + 1] - x) / h, b = (x - xs_[i]) / h;
    return a * ys_[i] + b * ys_[i + 1] + ((a * a * a - a) * m_[i] + (b * b * b - b) * m_[i + 1]) * (h * h / 6.0);
}

/** y'(x) = (y_{i+1} − y_i)/h − (3a² − 1) h M_i/6 + (3b² − 1) h M_{i+1}/6 */
template <Scalar T>
T NaturalCubic<T>::derivative(double x) const {
    if (extrap_ == Extrapolation::Flat && (x <= xs_.front() || x >= xs_.back())) return T(0.0);
    const double xc = std::min(std::max(x, xs_.front()), xs_.back()); // linear extrapolation continues the end slope
    const std::size_t i = detail::segment(xs_, xc);
    const double h = xs_[i + 1] - xs_[i];
    const double a = (xs_[i + 1] - xc) / h, b = (xc - xs_[i]) / h;
    return (ys_[i + 1] - ys_[i]) / h - (3.0 * a * a - 1.0) * h / 6.0 * m_[i] + (3.0 * b * b - 1.0) * h / 6.0 * m_[i + 1];
}

/** y''(x) = a M_i + b M_{i+1} (linear between nodes, zero at the ends: the "natural" condition). */
template <Scalar T>
T NaturalCubic<T>::secondDerivative(double x) const {
    if (x <= xs_.front() || x >= xs_.back()) return T(0.0);
    const std::size_t i = detail::segment(xs_, x);
    const double h = xs_[i + 1] - xs_[i];
    const double a = (xs_[i + 1] - x) / h, b = (x - xs_[i]) / h;
    return a * m_[i] + b * m_[i + 1];
}

} // namespace rl::math

// Explicit instantiation helper used by the two instantiation TUs.
#define RL_MATH_INSTANTIATE_INTERPOLATION(T)                                                           \
    template class rl::math::Linear<T>;                                                                \
    template class rl::math::LogLinearDf<T>;                                                           \
    template class rl::math::MonotoneConvex<T>;                                                        \
    template class rl::math::NaturalCubic<T>;                                                          \
    template T rl::math::interpolate<T>(const rl::math::Interpolator<T>&, double);                     \
    template T rl::math::interpolateDerivative<T>(const rl::math::Interpolator<T>&, double);
