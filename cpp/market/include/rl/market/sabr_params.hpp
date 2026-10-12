#pragma once
// SabrParams<T>: the parameters of one (shifted) SABR smile. A data record; the formula is vol/sabr.hpp.
//
//   dF̃ = α_t F̃^β dW,   dα = ν α dZ,   dW dZ = ρ dt,   F̃ = F + s        (Hagan et al. 2002, §3)
//
// Active set: α, ρ, ν are T — vega by cube node  differentiates through
// them. β and the shift s are passive doubles: they are desk *policy* chosen once per currency (β = 0.5 and
// s = 0 for USD in v1; β = 0 or a 1–3% shift where rates can go negative), never calibrated per node and
// never differentiated. Making them T would put policy on the tape and invite a fit to move them.
//
// Units: α has units of rate^{1−β} per √year (for β = 0.5 and F ≈ 3.5%, α ≈ 0.045); ρ is dimensionless;
// ν is per √year; s is a rate. Header-only: an aggregate and two inline functions.

#include "rl/core/error.hpp"
#include "rl/market/hash.hpp"
#include "rl/math/scalar.hpp"

#include <cmath>
#include <string>

namespace rl::market {

template <math::Scalar T>
struct SabrParams {
    T alpha = T(0.0);
    T rho = T(0.0);
    T nu = T(0.0);
    double beta = 0.5;
    double shift = 0.0;
};

/// Throws InputError unless α > 0, |ρ| < 1, ν ≥ 0, 0 ≤ β ≤ 1 and every field is finite.
template <math::Scalar T>
void validateSabrParams(const SabrParams<T>& p) {
    const double a = math::value(p.alpha), r = math::value(p.rho), n = math::value(p.nu);
    RL_REQUIRE(std::isfinite(a) && a > 0.0, "SABR alpha must be positive and finite, got " + std::to_string(a));
    RL_REQUIRE(std::isfinite(r) && r > -1.0 && r < 1.0, "SABR rho must lie in (-1, 1), got " + std::to_string(r));
    RL_REQUIRE(std::isfinite(n) && n >= 0.0, "SABR nu must be non-negative and finite, got " + std::to_string(n));
    RL_REQUIRE(std::isfinite(p.beta) && p.beta >= 0.0 && p.beta <= 1.0, "SABR beta must lie in [0, 1], got " + std::to_string(p.beta));
    RL_REQUIRE(std::isfinite(p.shift), "SABR shift must be finite");
}

/// Deterministic hash contribution (values, not tape nodes).
template <math::Scalar T>
void mixSabrParams(Hasher& h, const SabrParams<T>& p) noexcept {
    h.mix(math::value(p.alpha)).mix(math::value(p.rho)).mix(math::value(p.nu)).mix(p.beta).mix(p.shift);
}

} // namespace rl::market
