#pragma once
// definitions of the shifted Black-76 formulas (sources in black.hpp).
#include "rl/core/error.hpp"
#include "rl/math/normal.hpp"
#include "rl/vol/black.hpp"

#include <cmath>

namespace rl::vol {

namespace detail {
inline void requireBlackInputs(double shiftedForward, double expiry, double vol) {
    RL_REQUIRE(std::isfinite(shiftedForward) && shiftedForward > 0.0, "shifted Black: forward + shift must be positive");
    RL_REQUIRE(std::isfinite(expiry) && expiry >= 0.0, "shifted Black: expiry must be non-negative");
    RL_REQUIRE(std::isfinite(vol) && vol >= 0.0, "shifted Black: vol must be non-negative");
}
} // namespace detail

template <math::Scalar T>
T blackPrice(const T& forward, double strike, double expiry, const T& lognormalVol, double shift, OptionType type) {
    using std::log;
    const T f = forward + shift;
    const double k = strike + shift;
    detail::requireBlackInputs(math::value(f), expiry, math::value(lognormalVol));
    if (k <= 0.0) return type == OptionType::Call ? T(f - k) : T(0.0); // exercised almost surely / never
    const double w = optionSign(type);
    const T sd = lognormalVol * std::sqrt(expiry);
    if (math::value(sd) <= 0.0) return math::positivePart(T(w * (f - k)));
    const T d1 = log(f / k) / sd + 0.5 * sd;
    const T d2 = d1 - sd;
    return w * (f * math::normalCdf(T(w * d1)) - k * math::normalCdf(T(w * d2)));
}

template <math::Scalar T>
T blackDelta(const T& forward, double strike, double expiry, const T& lognormalVol, double shift, OptionType type) {
    using std::log;
    const T f = forward + shift;
    const double k = strike + shift;
    detail::requireBlackInputs(math::value(f), expiry, math::value(lognormalVol));
    if (k <= 0.0) return T(type == OptionType::Call ? 1.0 : 0.0);
    const double w = optionSign(type);
    const T sd = lognormalVol * std::sqrt(expiry);
    if (math::value(sd) <= 0.0) return T(w * (math::value(f) - k) > 0.0 ? w : 0.0);
    const T d1 = log(f / k) / sd + 0.5 * sd;
    return w * math::normalCdf(T(w * d1));
}

template <math::Scalar T>
T blackVega(const T& forward, double strike, double expiry, const T& lognormalVol, double shift) {
    using std::log;
    const T f = forward + shift;
    const double k = strike + shift;
    detail::requireBlackInputs(math::value(f), expiry, math::value(lognormalVol));
    if (k <= 0.0) return T(0.0);
    const double sqrtT = std::sqrt(expiry);
    const T sd = lognormalVol * sqrtT;
    if (math::value(sd) <= 0.0) return T(math::value(f) == k ? k * sqrtT * detail::kInvSqrt2Pi : 0.0);
    const T d1 = log(f / k) / sd + 0.5 * sd;
    return f * math::normalPdf(d1) * sqrtT;
}

template <math::Scalar T>
T blackGamma(const T& forward, double strike, double expiry, const T& lognormalVol, double shift) {
    using std::log;
    const T f = forward + shift;
    const double k = strike + shift;
    detail::requireBlackInputs(math::value(f), expiry, math::value(lognormalVol));
    if (k <= 0.0) return T(0.0);
    const T sd = lognormalVol * std::sqrt(expiry);
    if (math::value(sd) <= 0.0) return T(0.0);
    const T d1 = log(f / k) / sd + 0.5 * sd;
    return math::normalPdf(d1) / (f * sd);
}

} // namespace rl::vol

#define RL_VOL_INSTANTIATE_BLACK(T)                                                                             \
    namespace rl::vol {                                                                                         \
    template T blackPrice<T>(const T&, double, double, const T&, double, OptionType);                           \
    template T blackDelta<T>(const T&, double, double, const T&, double, OptionType);                           \
    template T blackVega<T>(const T&, double, double, const T&, double);                                        \
    template T blackGamma<T>(const T&, double, double, const T&, double);                                       \
    }
