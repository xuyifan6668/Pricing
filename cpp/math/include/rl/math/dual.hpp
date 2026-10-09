#pragma once
// Dual<T, N>: a forward-mode dual number with N tangent directions over a base scalar T ∈ {double, Real}.
// (Griewank & Walther, "Evaluating Derivatives" 2e, §3.1).
//
// What it is for: small, local Jacobians of templated formulas — the SABR smile's gradient in (α, ρ, ν)
//, the LGM strip's ζ sensitivities. Instantiating the *same* template on Dual<double, 3>
// gives the analytic gradient to rounding, with no hand-derived expression that can drift from the value code.
//
// DESIGN: D3 rejected forward-mode duals as the *risk engine* (O(n) passes for an n-node risk vector). That
// argument does not apply to a 3-parameter gradient evaluated at 7 strikes: forward mode is the cheapest mode
// there (one pass, no tape, no allocation), and it is exact. The tape stays the risk engine; Dual is a local
// tool. Dual<Real, N> composes: its value and tangents are Real, so a tape above it still records.
//
// DESIGN: satisfies rl::math::Scalar, and like rl::aad::Real it has NO relational operators: a branch
// inside a templated formula is written on value(x), so the tangent of a branch is the tangent of the branch
// taken — the almost-everywhere derivative, which is what both AAD and a Jacobian for a fit want.
//
// Header-only: a class template over two scalars and a handful of N; nothing to instantiate explicitly.

#include "rl/math/scalar.hpp"

#include <array>
#include <cstddef>

namespace rl::math {

template <Scalar T, std::size_t N>
class Dual {
public:
    using Base = T;
    static constexpr std::size_t kDirections = N;

    Dual() : v_(0.0) { d_.fill(T(0.0)); }
    Dual(double c) : v_(c) { d_.fill(T(0.0)); } // NOLINT: implicit by design (constants, as for aad::Real)
    Dual(const T& v, const std::array<T, N>& d) : v_(v), d_(d) {}

    /// A constant of the base type (zero tangent). For T = double this is the same as Dual(double).
    static Dual constant(const T& v) {
        std::array<T, N> d;
        d.fill(T(0.0));
        return Dual(v, d);
    }
    /// The i-th independent variable: value v, tangent e_i.
    static Dual variable(const T& v, std::size_t i) {
        Dual x = constant(v);
        x.d_[i] = T(1.0);
        return x;
    }

    const T& base() const noexcept { return v_; }
    const T& tangent(std::size_t i) const noexcept { return d_[i]; }
    const std::array<T, N>& tangents() const noexcept { return d_; }
    double value() const { return math::value(v_); }

    Dual& operator+=(const Dual& b) { return *this = *this + b; }
    Dual& operator-=(const Dual& b) { return *this = *this - b; }
    Dual& operator*=(const Dual& b) { return *this = *this * b; }
    Dual& operator/=(const Dual& b) { return *this = *this / b; }

    // --- arithmetic -----------------------------------------------------------------------------------
    friend Dual operator+(const Dual& a, const Dual& b) { return Dual(a.v_ + b.v_, zip(a.d_, b.d_, [](const T& x, const T& y) { return x + y; })); }
    friend Dual operator-(const Dual& a, const Dual& b) { return Dual(a.v_ - b.v_, zip(a.d_, b.d_, [](const T& x, const T& y) { return x - y; })); }
    friend Dual operator*(const Dual& a, const Dual& b) {
        return Dual(a.v_ * b.v_, zip(a.d_, b.d_, [&](const T& x, const T& y) { return x * b.v_ + a.v_ * y; }));
    }
    friend Dual operator/(const Dual& a, const Dual& b) {
        const T q = a.v_ / b.v_;
        return Dual(q, zip(a.d_, b.d_, [&](const T& x, const T& y) { return (x - q * y) / b.v_; }));
    }
    friend Dual operator+(const Dual& a, double b) { return Dual(a.v_ + b, a.d_); }
    friend Dual operator+(double a, const Dual& b) { return Dual(a + b.v_, b.d_); }
    friend Dual operator-(const Dual& a, double b) { return Dual(a.v_ - b, a.d_); }
    friend Dual operator-(double a, const Dual& b) { return Dual(a - b.v_, map(b.d_, [](const T& y) { return -y; })); }
    friend Dual operator*(const Dual& a, double b) { return Dual(a.v_ * b, map(a.d_, [&](const T& x) { return x * b; })); }
    friend Dual operator*(double a, const Dual& b) { return Dual(a * b.v_, map(b.d_, [&](const T& y) { return a * y; })); }
    friend Dual operator/(const Dual& a, double b) { return Dual(a.v_ / b, map(a.d_, [&](const T& x) { return x / b; })); }
    friend Dual operator/(double a, const Dual& b) {
        const T q = a / b.v_;
        return Dual(q, map(b.d_, [&](const T& y) { return -q * y / b.v_; }));
    }
    friend Dual operator-(const Dual& a) { return Dual(-a.v_, map(a.d_, [](const T& x) { return -x; })); }
    friend Dual operator+(const Dual& a) { return a; }

    // --- elementary functions (chain rule on the base scalar) --------------------------------------------
    friend Dual exp(const Dual& a) {
        const T e = exp(a.v_);
        return Dual(e, map(a.d_, [&](const T& x) { return e * x; }));
    }
    friend Dual log(const Dual& a) { return Dual(log(a.v_), map(a.d_, [&](const T& x) { return x / a.v_; })); }
    friend Dual sqrt(const Dual& a) {
        const T s = sqrt(a.v_);
        return Dual(s, map(a.d_, [&](const T& x) { return 0.5 * x / s; }));
    }
    /// d|a| = sign(a) da, with sign(0) = 0 (the subgradient convention of aad::Real).
    friend Dual fabs(const Dual& a) {
        const double sg = math::sign(a.v_);
        return Dual(fabs(a.v_), map(a.d_, [&](const T& x) { return sg * x; }));
    }
    friend Dual abs(const Dual& a) { return fabs(a); }
    friend Dual erf(const Dual& a) {
        const T g = kTwoOverSqrtPi * exp(-(a.v_ * a.v_));
        return Dual(erf(a.v_), map(a.d_, [&](const T& x) { return g * x; }));
    }
    friend Dual erfc(const Dual& a) {
        const T g = -kTwoOverSqrtPi * exp(-(a.v_ * a.v_));
        return Dual(erfc(a.v_), map(a.d_, [&](const T& x) { return g * x; }));
    }
    friend Dual pow(const Dual& a, double b) {
        const T g = b * pow(a.v_, b - 1.0);
        return Dual(pow(a.v_, b), map(a.d_, [&](const T& x) { return g * x; }));
    }
    /// ∂/∂b = a^b log a needs a > 0; for a ≤ 0 with an active exponent the tangent is NaN, deliberately.
    friend Dual pow(const Dual& a, const Dual& b) {
        const T p = pow(a.v_, b.v_);
        const T ga = b.v_ * pow(a.v_, b.v_ - 1.0);
        const T gb = p * log(a.v_);
        return Dual(p, zip(a.d_, b.d_, [&](const T& x, const T& y) { return ga * x + gb * y; }));
    }
    friend Dual pow(double a, const Dual& b) {
        const T p = pow(T(a), b.v_);
        const T gb = p * std::log(a);
        return Dual(p, map(b.d_, [&](const T& y) { return gb * y; }));
    }

private:
    static constexpr double kTwoOverSqrtPi = 1.12837916709551257390;

    template <class F>
    static std::array<T, N> map(const std::array<T, N>& x, F&& f) {
        std::array<T, N> out;
        for (std::size_t i = 0; i < N; ++i) out[i] = f(x[i]);
        return out;
    }
    template <class F>
    static std::array<T, N> zip(const std::array<T, N>& x, const std::array<T, N>& y, F&& f) {
        std::array<T, N> out;
        for (std::size_t i = 0; i < N; ++i) out[i] = f(x[i], y[i]);
        return out;
    }

    T v_;
    std::array<T, N> d_;
};

} // namespace rl::math
