#pragma once
// Bermudan swaption by finite differences.
//
// The deflated value V(ζ, x) is a martingale in a state with dx = sqrt(ζ') dW, so in the clock ζ
//
//     ∂V/∂ζ + ½ ∂²V/∂x² = 0,
//
// the heat equation with no drift, no discounting and constant coefficients: everything
// model-specific has been absorbed into the change of numeraire and the change of clock. The
// backward induction is
//     at the last exercise date:   V = max(h_n, 0)
//     between exercise dates:      roll V back with the heat equation
//     at each exercise date:       V ← max(V, h_i)
// and the price is V(0, 0).
//
// Scheme. Crank–Nicolson (second order in both step sizes), on a uniform grid over ±`stdDevs`
// standard deviations of the state at the last exercise date, with V_xx = 0 at the edges.
//
// Rannacher start-up, measured rather than assumed. Crank–Nicolson does not damp high-frequency
// error, and every max() puts a kink into V. The standard remedy replaces the first two steps after
// each exercise date by four fully implicit half-steps. On this product (docs/RESULTS.md):
//   * the PRICE is more accurate WITHOUT it at a given step count: the Bermudan's max() creates a
//     kink (a jump in the first derivative), not a jump in the value, which Crank–Nicolson handles
//     well, and the implicit steps add their own first-order-in-time error;
//   * GAMMA is 3–6 times smoother WITH it: the undamped modes show up in the second derivative.
// It is on by default because a pricer's output is hedged, not just quoted.
//
// Space error is second order on average but not monotone in the step: it depends on where the
// exercise boundary falls between grid nodes, so refining can leave the error unchanged before the
// next refinement cuts it fifty-fold. Convergence is therefore judged over several doublings.

#include "rl/pricers/bermudan.hpp"

#include <vector>

namespace rl::pricers {

struct PdeOptions {
    int spaceSteps = 400;    ///< grid has spaceSteps + 1 points (made even so that x = 0 is a node)
    int timeSteps = 200;     ///< total over [0, ζ_n], split across intervals in proportion to Δζ (≥ 4 each)
    double stdDevs = 7.5;
    bool rannacher = true;
};

struct PdeResult {
    double value = 0.0;
    /// State at which exercise becomes optimal at each date (NaN if exercise is never / always optimal on the grid).
    std::vector<double> boundaryState;
    /// The same boundary expressed as the par rate of the swap entered into.
    std::vector<double> boundaryRate;
    std::vector<double> grid;   ///< x
    std::vector<double> values; ///< V(0, x): today's value as a function of a hypothetical starting state
    int spaceSteps = 0;
    int timeSteps = 0;
};

PdeResult priceBermudanPde(const BermudanProblem& problem, const PdeOptions& options = {});

} // namespace rl::pricers
