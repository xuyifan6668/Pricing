#pragma once
// Does the delta hedge work? A hedging simulation and a P&L explain for the Bermudan.
//
// HEDGING SIMULATION. A desk sells the Bermudan and delta-hedges it with the swap it would deliver
// on the first exercise date, rebalancing at a fixed frequency until that date, when the position
// is closed at its model value. In a one-factor model with one hedge instrument the market is
// complete: with continuous rebalancing and the right model the hedge replicates the option and the
// P&L is zero on every path. The simulation measures how far from zero it is
//
//   * with discrete rebalancing: the error is gamma times the squared move over each interval, so
//     its standard deviation falls like the square root of the interval;
//   * with the wrong volatility: the hedger computes its price and deltas from a model whose
//     variance is scaled by s², while the market moves with the true one. The mean P&L is then the
//     mispricing, and the error no longer vanishes as rebalancing gets faster.
//
// Everything is in deflated units, where both the option and the hedge swap are martingales, so
// there is no financing term: the hedge gain over an interval is Δ times the change in the swap's
// deflated value. Option values between exercise dates come from the value on the first exercise
// date's grid integrated against the exact Gaussian transition density; deltas from the same
// integral differentiated in the state.
//
// The simulation stops at the first exercise date. Hedging through the exercise decisions would add
// the question of what the hedger does when the counterparty exercises, which this does not model.

#include "rl/models/usd_market.hpp"
#include "rl/pricers/bermudan.hpp"

#include <cstddef>
#include <cstdint>
#include <functional>
#include <vector>

namespace rl::pricers {

struct HedgingOptions {
    std::size_t paths = 20000;
    int rebalancesPerYear = 52;
    double volScale = 1.0;      ///< the hedger's model has variance volScale² times the true one
    bool hedge = true;          ///< false: sell and hold, for comparison
    std::uint64_t seed = 20260904;
    int gridPoints = 801;
};

struct HedgingResult {
    double premium = 0.0;        ///< the hedger's price today
    double truePrice = 0.0;      ///< the true model's price today
    double meanPnl = 0.0;        ///< of the short-and-hedged position, closed at the first exercise date
    double pnlStdDev = 0.0;
    double meanPnlStandardError = 0.0;
    double horizon = 0.0;        ///< years: the first exercise date
};

/// `trueModel` drives the market; the hedger uses the same mean reversion with variance scaled by
/// `volScale²`. `rights` are the Bermudan's exercise rights.
HedgingResult simulateDeltaHedge(const models::Lgm1F& trueModel, const std::vector<models::ExerciseRight>& rights, const HedgingOptions& options = {});

// --------------------------------------------------------------------------------------------------
// P&L explain
// --------------------------------------------------------------------------------------------------

/// A market move: a shift in bp for each curve quote (rates up; a future's price down by 0.01 per
/// bp, as in risk.hpp) and a parallel shift of every swaption vol quote, in bp of normal vol.
struct MarketMove {
    std::vector<double> quoteShiftBp;
    double volShiftBp = 0.0;
};

/// The market after a move: curve quotes shifted and re-bootstrapped, vol screen shifted and refitted.
models::UsdMarket applyMove(const models::UsdMarket& market, const MarketMove& move);

/// The sensitivities a risk system would report for a price: par-quote delta (per bp, central ±1 bp
/// bumps of each quote), parallel gamma (per bp², from ±`gammaBumpBp` parallel shifts) and parallel
/// vega (per bp of normal vol, ±1 bp central).
///
/// Why gamma takes a wider bump. A second difference divides by the bump squared, so any noise in the
/// price is amplified by 1/h². A PDE price re-solved after re-calibration moves its grid with the
/// model's variance, which leaves a numerical wiggle of a few thousandths of a basis point: on a 400-step
/// grid a ±1 bp second difference of the Bermudan is pure noise (it even has the wrong sign), while
/// ±10 bp recovers the convexity. The same reasoning is why production systems fix bump sizes per risk
/// type and test them, rather than using the smallest bump that runs.
struct Sensitivities {
    double value = 0.0;
    std::vector<double> delta;
    double parallelGamma = 0.0;
    double vega = 0.0;
};
Sensitivities computeSensitivities(const models::UsdMarket& market, const std::function<double(const models::UsdMarket&)>& price, double gammaBumpBp = 10.0);

struct PnlExplain {
    double actual = 0.0;      ///< full revaluation
    double delta = 0.0;       ///< Σ δ_i shift_i
    double gamma = 0.0;       ///< ½ Γ (mean shift)²: parallel gamma applied to the move's average
    double vega = 0.0;        ///< vega × vol shift
    double unexplained() const noexcept { return actual - delta - gamma - vega; }
};

/// Attribute the change in `price` under `move` to the sensitivities.
PnlExplain explainPnl(const models::UsdMarket& market, const Sensitivities& s, const std::function<double(const models::UsdMarket&)>& price, const MarketMove& move);

} // namespace rl::pricers
