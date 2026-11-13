#pragma once
// The SABR MODEL by simulation, against Hagan's FORMULA for it.
//
// Everything else in vol/ uses Hagan's normal-vol expansion: an asymptotic approximation to the
// model's implied vol, accurate when ν²T is small. The cube is fitted with the same formula, so its
// fit error says nothing about the expansion's error. This prices the model itself:
//
//     dF̃ = α F̃^β dW,     dα = ν α dZ,     dW dZ = ρ dt,     F̃ = F + shift,
//
// under the forward (annuity) measure, where F is a martingale and an option is E[(F_T − K)⁺] per
// unit of annuity. Simulation:
//   * α is simulated exactly: α_{t+Δ} = α_t exp(ν √Δ Z − ½ ν² Δ).
//   * F̃ by Euler with the α at the start of the step, absorbed at zero when β > 0 (the standard
//     boundary condition: a rate that reaches zero stays there). With β = 0 F is not floored.
//   * Antithetic pairs, and the terminal forward as a control variate (its mean is known: F̃_0). The
//     control also absorbs most of the drift Euler with absorption introduces.
//   * Out-of-the-money options only: puts below the forward, calls above, so that the implied vol is
//     read from time value, not from intrinsic value plus noise.
//
// What it is for. The true model has a non-negative density by construction; Hagan's formula does
// not (vol/arbitrage.hpp). Comparing the two node by node says where the expansion can be trusted,
// which the fit residual cannot.
//
// The Euler error is first order in the step: the tests and RESULTS compare step sizes.

#include "rl/market/sabr_params.hpp"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace rl::vol {

struct SabrMcOptions {
    std::size_t paths = 200000; ///< even: antithetic pairs
    int stepsPerYear = 100;
    std::uint64_t seed = 20260904;
};

struct SabrMcPoint {
    double strike = 0.0;
    double price = 0.0;          ///< out-of-the-money option per unit annuity (put below the forward, call at or above)
    double standardError = 0.0;
    double normalVol = 0.0;      ///< implied Bachelier vol of `price`; NaN if the price is too small to invert
    double normalVolError = 0.0; ///< standard error mapped to vol through vega
    double haganNormalVol = 0.0; ///< Hagan's expansion at the same strike
};

struct SabrMcResult {
    std::vector<SabrMcPoint> points;
    double forwardError = 0.0;   ///< E[F_T] − F_0 on the paths, before the control variate
    double absorbed = 0.0;       ///< fraction of paths that hit zero
};

/// Throws InputError on invalid parameters, a non-positive expiry, or an odd or small path count.
SabrMcResult sabrMonteCarlo(const market::SabrParams<double>& params, double forward, double expiry, const std::vector<double>& strikes, const SabrMcOptions& options = {});

} // namespace rl::vol
