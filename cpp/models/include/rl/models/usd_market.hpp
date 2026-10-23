#pragma once
// The market every example, test and analysis in this repository prices against: the USD curves
// bootstrapped from data/quotes, and the swaption cube fitted from data/vol, as of 2026-09-04.
// Building it runs the two calibrations (global bootstrap, 81 SABR fits); it takes a few milliseconds.

#include "rl/bootstrap/bootstrap_problem.hpp"
#include "rl/bootstrap/global_bootstrap.hpp"
#include "rl/market/curve_interpolation.hpp"
#include "rl/market/curve_set.hpp"
#include "rl/market/fixings.hpp"
#include "rl/market/swaption_cube.hpp"
#include "rl/market/vol_quotes.hpp"
#include "rl/vol/sabr_fit.hpp"

namespace rl::models {

struct UsdMarket {
    core::Date asOf;
    bootstrap::BootstrapProblem problem; ///< the curve quotes and instruments (for par-quote bumps)
    math::Vector z;                      ///< solved curve nodes
    market::CurveSet<double> curves;
    market::FixingStore fixings;
    market::SwaptionVolQuotes volQuotes; ///< the vol screen the cube was fitted to
    market::SwaptionCube<double> cube;
    vol::CubeFit cubeFit;                ///< per-node SABR fit results

    /// The same market with the curve quotes replaced (decimal rates; futures in price) and the
    /// curves re-bootstrapped. The vol cube is unchanged: this is a rates-only bump.
    UsdMarket withCurveQuotes(const math::Vector& quotes) const;
    /// The same market with every ATM normal vol and smile quote shifted by `bump` (decimal) and the
    /// cube re-fitted.
    UsdMarket withVolShift(double bump) const;
};

UsdMarket loadUsdMarket(market::CurveInterpolation interpolation = market::CurveInterpolation::LinearLogDf);

} // namespace rl::models
