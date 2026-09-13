#pragma once
// Interpolation strategies (docs/DESIGN.md §2.5: a std::variant, dispatched once per query;
// never inside a per-path loop, because the model precomputes on the timeline).
//
// Every interpolator maps a passive abscissa x (a time in years, ACT/365F) to an active ordinate
// y ∈ T. Node abscissae are strictly increasing doubles; node ordinates are T.
//
// Which quantity to interpolate is the curve's business. The intended pairing is:
//   Linear on y = −log P(0,t)     → piecewise-flat instantaneous forwards, the market default for
//                                   OIS curves: local, cheap, never oscillates, forwards jump at nodes.
//   LogLinearDf on y = P(0,t)      → the same thing expressed on discount factors (for the notebook).
//   MonotoneConvex on y = −log P   → Hagan–West 2006: continuous forwards, positivity preserved,
//                                   local (a bump moves at most two intervals). Production choice.
//   NaturalCubic on y = σ(t)       → the vol cube's expiry/tenor axes; NOT for curves (global support,
//                                   overshoots, forwards go negative).
//
// Thread-safety: immutable after construction; share freely.

#include "rl/math/scalar.hpp"

#include <cstddef>
#include <variant>
#include <vector>

namespace rl::math {

enum class Extrapolation {
    Flat,   ///< y(x) = y(x_0) for x < x_0, y(x_{n−1}) beyond.
    Linear, ///< Continue the first/last segment's slope (flat *forward* when y = −log P).
};

namespace detail {
/// Locate the segment i with x_i ≤ x < x_{i+1}; clamps to [0, n−2] outside the range. O(log n).
std::size_t segment(const std::vector<double>& xs, double x) noexcept;
/// Throws InputError unless xs is strictly increasing with ≥ 2 points and ys matches in size.
void validateNodes(const std::vector<double>& xs, std::size_t ny, std::size_t minPoints);
} // namespace detail

template <Scalar T>
class Linear {
public:
    Linear(std::vector<double> xs, std::vector<T> ys, Extrapolation extrap = Extrapolation::Linear);

    T operator()(double x) const;   ///< y(x)
    T derivative(double x) const;   ///< dy/dx (piecewise constant; at a node, the right-hand slope)

    const std::vector<double>& xs() const noexcept { return xs_; }
    const std::vector<T>& ys() const noexcept { return ys_; }
    std::size_t size() const noexcept { return xs_.size(); }

private:
    std::vector<double> xs_;
    std::vector<T> ys_;
    Extrapolation extrap_;
};

/// Interpolates log(y) linearly; ys must be strictly positive (discount factors).
template <Scalar T>
class LogLinearDf {
public:
    LogLinearDf(std::vector<double> xs, std::vector<T> dfs, Extrapolation extrap = Extrapolation::Linear);

    T operator()(double x) const;   ///< P(x)
    T derivative(double x) const;   ///< dP/dx = −f(x)·P(x)
    T forward(double x) const;      ///< instantaneous forward −d log P/dx (piecewise constant)

    const std::vector<double>& xs() const noexcept { return logInterp_.xs(); }
    std::size_t size() const noexcept { return logInterp_.size(); }

private:
    Linear<T> logInterp_; // on log(df)
};

/// Hagan & West (2006), "Interpolation Methods for Curve Construction", Applied Mathematical
/// Finance 13(2), §5–6 ("monotone convex"). Nodes are (t_i, y_i) with y_i = r_i·t_i = −log P(0,t_i),
/// t_0 = 0, y_0 = 0. The method interpolates the *instantaneous forward* f(t) such that
/// ∫_{t_{i−1}}^{t_i} f = y_i − y_{i−1}, preserving monotonicity and (where the data allow) positivity.
template <Scalar T>
class MonotoneConvex {
public:
    MonotoneConvex(std::vector<double> ts, std::vector<T> ys, bool enforcePositivity = true);

    T operator()(double t) const;   ///< y(t) = ∫_0^t f(s) ds
    T derivative(double t) const { return forward(t); }
    T forward(double t) const;      ///< f(t), the interpolated instantaneous forward

    const std::vector<double>& xs() const noexcept { return ts_; }
    const std::vector<T>& ys() const noexcept { return ys_; }
    /// Discrete forwards f^d_i on (t_{i−1}, t_i], i = 1..n−1 (size n−1).
    const std::vector<T>& discreteForwards() const noexcept { return fd_; }
    /// Node forwards f_i at t_i, i = 0..n−1 (size n), after the Hagan–West collar.
    const std::vector<T>& nodeForwards() const noexcept { return f_; }
    std::size_t size() const noexcept { return ts_.size(); }

private:
    void build(); // computes fd_ and f_ (the maths; see interpolation.ipp)
    std::vector<double> ts_;
    std::vector<T> ys_;
    std::vector<T> fd_;
    std::vector<T> f_;
    bool positivity_;
};

/// Natural cubic spline (second derivative zero at both ends). Global: every ordinate influences
/// every segment. Flat extrapolation by default (a vol surface should not be extrapolated with a cubic).
template <Scalar T>
class NaturalCubic {
public:
    NaturalCubic(std::vector<double> xs, std::vector<T> ys, Extrapolation extrap = Extrapolation::Flat);

    T operator()(double x) const;
    T derivative(double x) const;
    T secondDerivative(double x) const;

    const std::vector<double>& xs() const noexcept { return xs_; }
    const std::vector<T>& ys() const noexcept { return ys_; }
    /// Second derivatives at the nodes (M_i); M_0 = M_{n−1} = 0.
    const std::vector<T>& secondDerivatives() const noexcept { return m_; }
    std::size_t size() const noexcept { return xs_.size(); }

private:
    void build(); // solves the tridiagonal system for m_ (the maths)
    std::vector<double> xs_;
    std::vector<T> ys_;
    std::vector<T> m_;
    Extrapolation extrap_;
};

/// Runtime-selected strategy. Dispatched with std::visit once per query.
template <Scalar T>
using Interpolator = std::variant<Linear<T>, LogLinearDf<T>, MonotoneConvex<T>, NaturalCubic<T>>;

template <Scalar T>
T interpolate(const Interpolator<T>& interp, double x) {
    return std::visit([x](const auto& f) { return f(x); }, interp);
}

template <Scalar T>
T interpolateDerivative(const Interpolator<T>& interp, double x) {
    return std::visit([x](const auto& f) { return f.derivative(x); }, interp);
}

} // namespace rl::math
