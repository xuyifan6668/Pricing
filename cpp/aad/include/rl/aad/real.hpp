#pragma once
// rl::aad::Real — the active scalar. Satisfies rl::math::Scalar (see scalar.hpp) and records every
// operation on the thread-local Tape.
//
// Two kinds of Real:
//   constant  (node == Tape::kPassive): built from a double, records nothing. `Real x = 2.0;`
//   active    (node < tape size): a leaf made with Real::input(v), or the result of an operation
//             with at least one active operand.
// Operations between two constants produce a constant — so a pure-double computation wrapped in
// Real costs an extra branch per op and no tape traffic, which is the "identical values, no
// adjoint" mode the tests exercise.
//
// There are NO relational operators: write `value(x) < 0` (see scalar.hpp for why).
//
// Every operator records its local partial derivatives at evaluation time, so the reverse sweep
// (Tape::propagate) is a pure multiply-add over the tape with no knowledge of what the nodes were.

#include "rl/aad/tape.hpp"

#include <cmath>
#include <cstdint>
#include <iosfwd>
#include <limits>

namespace rl::aad {

inline constexpr double kTwoOverSqrtPi = 1.1283791670955126; // 2/√π, the derivative of erf at 0

class Real {
public:
    constexpr Real() noexcept : v_(0.0), node_(Tape::kPassive) {}
    constexpr Real(double v) noexcept : v_(v), node_(Tape::kPassive) {} // NOLINT: implicit by design (constants)

    /// An active leaf: registers an input on the current tape.
    static Real input(double v) { return Real(v, Tape::current().newLeaf()); }

    constexpr double value() const noexcept { return v_; }
    constexpr std::uint32_t node() const noexcept { return node_; }
    constexpr bool isActive() const noexcept { return node_ != Tape::kPassive; }

    // --- compound assignment --------------------------------------------------------------------
    Real& operator+=(const Real& b) { return *this = *this + b; }
    Real& operator-=(const Real& b) { return *this = *this - b; }
    Real& operator*=(const Real& b) { return *this = *this * b; }
    Real& operator/=(const Real& b) { return *this = *this / b; }
    Real& operator+=(double b) { return *this = *this + b; }
    Real& operator-=(double b) { return *this = *this - b; }
    Real& operator*=(double b) { return *this = *this * b; }
    Real& operator/=(double b) { return *this = *this / b; }

    // --- arithmetic; each records one node with the parents' partials ---
    /// ∂/∂a = 1, ∂/∂b = 1
    friend Real operator+(const Real& a, const Real& b) { return binary(a.v_ + b.v_, a, 1.0, b, 1.0); }
    /// ∂/∂a = 1, ∂/∂b = −1
    friend Real operator-(const Real& a, const Real& b) { return binary(a.v_ - b.v_, a, 1.0, b, -1.0); }
    /// ∂/∂a = b, ∂/∂b = a
    friend Real operator*(const Real& a, const Real& b) { return binary(a.v_ * b.v_, a, b.v_, b, a.v_); }
    /// ∂/∂a = 1/b, ∂/∂b = −a/b²
    friend Real operator/(const Real& a, const Real& b) { return binary(a.v_ / b.v_, a, 1.0 / b.v_, b, -a.v_ / (b.v_ * b.v_)); }

    /// ∂/∂a = 1
    friend Real operator+(const Real& a, double b) { return unary(a.v_ + b, a, 1.0); }
    friend Real operator+(double a, const Real& b) { return unary(a + b.v_, b, 1.0); }
    /// ∂/∂a = 1;  for (double − Real): ∂/∂b = −1
    friend Real operator-(const Real& a, double b) { return unary(a.v_ - b, a, 1.0); }
    friend Real operator-(double a, const Real& b) { return unary(a - b.v_, b, -1.0); }
    /// ∂/∂a = b (the constant)
    friend Real operator*(const Real& a, double b) { return unary(a.v_ * b, a, b); }
    friend Real operator*(double a, const Real& b) { return unary(a * b.v_, b, a); }
    /// ∂/∂a = 1/b;  for (double / Real): ∂/∂b = −a/b²
    friend Real operator/(const Real& a, double b) { return unary(a.v_ / b, a, 1.0 / b); }
    friend Real operator/(double a, const Real& b) { return unary(a / b.v_, b, -a / (b.v_ * b.v_)); }

    /// ∂/∂a = −1
    friend Real operator-(const Real& a) { return unary(-a.v_, a, -1.0); }
    friend Real operator+(const Real& a) { return a; }

    // --- elementary functions ------------------------------------------------------------------
    /// ∂/∂a = exp(a)   (reuse the result value)
    friend Real exp(const Real& a) {
        const double e = std::exp(a.v_);
        return unary(e, a, e);
    }
    /// ∂/∂a = 1/a
    friend Real log(const Real& a) { return unary(std::log(a.v_), a, 1.0 / a.v_); }
    /// ∂/∂a = 1/(2 √a)   (reuse the result; at a = 0 the derivative is +∞ — document, do not guard)
    friend Real sqrt(const Real& a) {
        const double r = std::sqrt(a.v_);
        return unary(r, a, 0.5 / r);
    }
    /// ∂/∂a = sign(a)   (0 at a = 0 by convention)
    friend Real fabs(const Real& a) { return unary(std::fabs(a.v_), a, static_cast<double>((a.v_ > 0.0) - (a.v_ < 0.0))); }
    friend Real abs(const Real& a) { return fabs(a); }
    /// ∂/∂a = (2/√π) exp(−a²)
    friend Real erf(const Real& a) { return unary(std::erf(a.v_), a, kTwoOverSqrtPi * std::exp(-a.v_ * a.v_)); }
    /// ∂/∂a = −(2/√π) exp(−a²)
    friend Real erfc(const Real& a) { return unary(std::erfc(a.v_), a, -kTwoOverSqrtPi * std::exp(-a.v_ * a.v_)); }
    /// ∂/∂a = b a^{b−1},  ∂/∂b = a^b log a   (log a needs a > 0; for a ≤ 0 with active b the
    /// derivative is undefined — record NaN deliberately, that is a modelling error upstream)
    friend Real pow(const Real& a, const Real& b) {
        const double p = std::pow(a.v_, b.v_);
        return binary(p, a, b.v_ * std::pow(a.v_, b.v_ - 1.0), b, p * std::log(a.v_));
    }
    /// ∂/∂a = b a^{b−1}
    friend Real pow(const Real& a, double b) { return unary(std::pow(a.v_, b), a, b * std::pow(a.v_, b - 1.0)); }
    /// ∂/∂b = a^b log a
    friend Real pow(double a, const Real& b) {
        const double p = std::pow(a, b.v_);
        return unary(p, b, p * std::log(a));
    }

    /// Equality on values only. Eigen's product kernels need == on the scalar; ordering stays absent
    /// so a branch on an active variable still has to be written on value(x).
    friend constexpr bool operator==(const Real& a, const Real& b) noexcept { return a.v_ == b.v_; }
    friend constexpr bool operator!=(const Real& a, const Real& b) noexcept { return a.v_ != b.v_; }

    /// Passive value (the Scalar concept's value()).
    friend constexpr double value(const Real& a) noexcept { return a.v_; }

private:
    constexpr Real(double v, std::uint32_t node) noexcept : v_(v), node_(node) {}

    /// Record a one-parent node unless the parent is a constant (then the result is a constant).
    static Real unary(double v, const Real& a, double da) {
        if (!a.isActive()) return Real(v);
        return Real(v, Tape::current().record1(a.node_, da));
    }
    /// Record a two-parent node; constants among the parents are dropped by the tape.
    static Real binary(double v, const Real& a, double da, const Real& b, double db) {
        if (!a.isActive() && !b.isActive()) return Real(v);
        return Real(v, Tape::current().record2(a.node_, da, b.node_, db));
    }

    double v_;
    std::uint32_t node_;
};

std::ostream& operator<<(std::ostream& os, const Real& r);

} // namespace rl::aad
