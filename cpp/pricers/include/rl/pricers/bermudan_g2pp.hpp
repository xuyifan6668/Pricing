#pragma once
// The Bermudan swaption under the two-factor model, by least-squares Monte Carlo.
//
// With two state variables a finite-difference grid is two-dimensional and Monte Carlo becomes the
// natural method. The estimator is the Longstaff–Schwartz lower bound of bermudan_mc.hpp with the
// two changes the second factor forces:
//
//   * The state is simulated exactly. Under the numeraire P(·, T*), T* the last cashflow date, the
//     state has no drift and its increment between exercise dates is N(0, C(t_i) − C(t_{i−1})), so
//     there is no time-stepping error. T* is chosen because under the model's textbook numeraire a
//     fast factor makes path weights vary over many orders of magnitude (g2pp.hpp).
//
//   * The continuation value is regressed on polynomials in BOTH states. With one factor the swap's
//     value is a monotone function of the state, so "the state" and "how far in the money" are the
//     same information. With two, paths with the same exercise value can have different curve
//     shapes, and so different continuation values; a regression on the exercise value alone would
//     throw that away. The regressors are the first state and the part of the second that is
//     independent of it, both standardised.
//
// What this is for. Calibrated to the same co-terminal strip, the one- and two-factor models agree
// on every European the Bermudan contains. They disagree on how correlated the swap rates seen on
// different exercise dates are, and that is exactly what the switch value prices. Running this
// pricer with the second factor switched off reproduces the one-factor price, so the difference
// between the two runs is the effect of the second factor and nothing else.

#include "rl/models/flows.hpp"
#include "rl/models/g2pp.hpp"
#include "rl/models/g2pp_calibration.hpp"
#include "rl/models/usd_market.hpp"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace rl::pricers {

struct G2ppLsmcOptions {
    std::size_t trainingPaths = 20000;
    std::size_t pricingPaths = 100000; ///< independent of the training paths
    int degree = 3;                    ///< total degree of the polynomial in the two states
    bool antithetic = true;            ///< pricing pass only
    std::uint64_t seed = 20260904;
};

struct G2ppLsmcResult {
    double lower = 0.0;         ///< out-of-sample: a valid lower bound up to its standard error
    double standardError = 0.0;
    double inSample = 0.0;      ///< the rule on its own training paths: biased high
    std::vector<double> exerciseProbability; ///< fraction of pricing paths exercised at each date
    std::vector<double> europeans;           ///< model price of the European with each single exercise date
};

/// Rights must be in increasing exercise time. Throws InputError otherwise.
G2ppLsmcResult priceBermudanLsmcG2pp(const models::G2pp& model, const std::vector<models::ExerciseRight>& rights, const G2ppLsmcOptions& options = {});

/// A Bermudan with a two-factor model whose volatility scale is calibrated to its own co-terminal strip.
struct CalibratedBermudanG2pp {
    std::vector<models::ExerciseRight> rights;
    models::G2ppStripCalibration calibration;
};

/// Exercise rights from `swap` (annual call dates after `firstCallYears`, `noticeDays` business days
/// of notice), strip from the market's cube at the swap's fixed rate, scale of `shape` bootstrapped to it.
CalibratedBermudanG2pp calibrateBermudanG2pp(const models::UsdMarket& market, const instruments::Swap& swap, double firstCallYears, int noticeDays,
                                             const models::G2pp& shape);

} // namespace rl::pricers
