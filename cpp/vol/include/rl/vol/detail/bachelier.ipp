#pragma once
// definitions of the Bachelier formulas (sources in bachelier.hpp). Included by src/vol.cpp (double)
// and src/instantiate_real.cpp (Real).
#include "rl/core/error.hpp"
#include "rl/math/normal.hpp"
#include "rl/vol/bachelier.hpp"

#include <cmath>

namespace rl::vol {

namespace detail {
inline constexpr double kInvSqrt2Pi = 0.39894228040143267794;

inline void requireBachelierInputs(double expiry, double vol) {
    RL_REQUIRE(std::isfinite(expiry) && expiry >= 0.0, "Bachelier: expiry must be non-negative");
    RL_REQUIRE(std::isfinite(vol) && vol >= 0.0, "Bachelier: normal vol must be non-negative");
}
} // namespace detail

template <math::Scalar T>
T bachelierPrice(const T& forward, double strike, double expiry, const T& normalVol, OptionType type) {
    detail::requireBachelierInputs(expiry, math::value(normalVol));
    const double w = optionSign(type);
    const T m = forward - strike;
    const T sd = normalVol * std::sqrt(expiry);
    if (math::value(sd) <= 0.0) return math::positivePart(T(w * m)); // no time value
    const T d = m / sd;
    return w * m * math::normalCdf(T(w * d)) + sd * math::normalPdf(d);
}

template <math::Scalar T>
T bachelierDelta(const T& forward, double strike, double expiry, const T& normalVol, OptionType type) {
    detail::requireBachelierInputs(expiry, math::value(normalVol));
    const double w = optionSign(type);
    const T m = forward - strike;
    const T sd = normalVol * std::sqrt(expiry);
    if (math::value(sd) <= 0.0) return T(w * math::value(m) > 0.0 ? w : 0.0);
    return w * math::normalCdf(T(w * (m / sd)));
}

template <math::Scalar T>
T bachelierVega(const T& forward, double strike, double expiry, const T& normalVol) {
    detail::requireBachelierInputs(expiry, math::value(normalVol));
    const double sqrtT = std::sqrt(expiry);
    const T m = forward - strike;
    const T sd = normalVol * sqrtT;
    if (math::value(sd) <= 0.0) return T(math::value(m) == 0.0 ? sqrtT * detail::kInvSqrt2Pi : 0.0);
    return sqrtT * math::normalPdf(T(m / sd));
}

template <math::Scalar T>
T bachelierGamma(const T& forward, double strike, double expiry, const T& normalVol) {
    detail::requireBachelierInputs(expiry, math::value(normalVol));
    const T sd = normalVol * std::sqrt(expiry);
    if (math::value(sd) <= 0.0) return T(0.0);
    return math::normalPdf(T((forward - strike) / sd)) / sd;
}

} // namespace rl::vol

// One line per template; a missing line is an unresolved external.
#define RL_VOL_INSTANTIATE_BACHELIER(T)                                                                         \
    namespace rl::vol {                                                                                         \
    template T bachelierPrice<T>(const T&, double, double, const T&, OptionType);                               \
    template T bachelierDelta<T>(const T&, double, double, const T&, OptionType);                               \
    template T bachelierVega<T>(const T&, double, double, const T&);                                            \
    template T bachelierGamma<T>(const T&, double, double, const T&);                                           \
    }
