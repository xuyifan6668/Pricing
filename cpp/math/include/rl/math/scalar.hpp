#pragma once
// The Scalar concept: what every templated numerical routine may assume about T.
//
// T ∈ {double, rl::aad::Real}. Templates downstream write `using std::exp; exp(x)` (or call through
// this namespace, which does that for them) so that ADL finds the aad overloads for Real and the
// std ones for double.
//
// DELIBERATE OMISSION: relational operators are not part of the concept and rl::aad::Real does not
// define them. A branch on an active variable is written `if (value(x) < 0)`, which makes the
// active-set decision visible at every branch and turns the classic
// AAD bug — control flow silently baked at record time — into a compile error instead. max/min are
// provided here as free functions that compare values and return one of the inputs (so the tape sees
// an identity, which is the correct almost-everywhere derivative).

#include <cmath>
#include <concepts>
#include <type_traits>

namespace rl::math {

using std::erf;
using std::erfc;
using std::exp;
using std::fabs;
using std::log;
using std::pow;
using std::sqrt;

/// Passive value of a scalar: identity on arithmetic types, x.value() on an active type. Works
/// both unqualified (ADL also finds rl::aad::value) and as rl::math::value(x).
template <class T>
constexpr double value(const T& x) noexcept {
    if constexpr (std::is_arithmetic_v<T>) return static_cast<double>(x);
    else return x.value();
}

template <class T>
concept Scalar = requires(T a, T b, double d) {
    T{};
    T(d);
    { a + b } -> std::convertible_to<T>;
    { a - b } -> std::convertible_to<T>;
    { a * b } -> std::convertible_to<T>;
    { a / b } -> std::convertible_to<T>;
    { a + d } -> std::convertible_to<T>;
    { d + a } -> std::convertible_to<T>;
    { a - d } -> std::convertible_to<T>;
    { d - a } -> std::convertible_to<T>;
    { a * d } -> std::convertible_to<T>;
    { d * a } -> std::convertible_to<T>;
    { a / d } -> std::convertible_to<T>;
    { d / a } -> std::convertible_to<T>;
    { -a } -> std::convertible_to<T>;
    { a += b } -> std::same_as<T&>;
    { a -= b } -> std::same_as<T&>;
    { a *= b } -> std::same_as<T&>;
    { a /= b } -> std::same_as<T&>;
    { exp(a) } -> std::convertible_to<T>;
    { log(a) } -> std::convertible_to<T>;
    { sqrt(a) } -> std::convertible_to<T>;
    { fabs(a) } -> std::convertible_to<T>;
    { erf(a) } -> std::convertible_to<T>;
    { erfc(a) } -> std::convertible_to<T>;
    { pow(a, b) } -> std::convertible_to<T>;
    { pow(a, d) } -> std::convertible_to<T>;
    { value(a) } -> std::convertible_to<double>;
};

static_assert(Scalar<double>);

/// max/min by value comparison, returning one of the inputs (see header comment).
template <Scalar T>
constexpr const T& max(const T& a, const T& b) noexcept { return value(a) < value(b) ? b : a; }
template <Scalar T>
constexpr const T& min(const T& a, const T& b) noexcept { return value(b) < value(a) ? b : a; }

/// (x)^+ = max(x, 0) with the zero as a T so the return type is T.
template <Scalar T>
T positivePart(const T& x) { return value(x) > 0.0 ? x : T(0.0); }

template <Scalar T>
bool isFinite(const T& x) noexcept { return std::isfinite(value(x)); }

/// sign(x) ∈ {−1, 0, +1} as a double (passive by nature).
template <Scalar T>
constexpr double sign(const T& x) noexcept {
    const double v = value(x);
    return (v > 0.0) - (v < 0.0);
}

} // namespace rl::math
