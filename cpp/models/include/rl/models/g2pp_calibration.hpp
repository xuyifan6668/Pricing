#pragma once
// Calibrating G2++ in two steps, which answer two different questions.
//
// 1. SHAPE, from the at-the-money swaption matrix. The matrix has two directions. Along EXPIRY the
//    level of volatility changes, which any Gaussian model absorbs with a time-dependent volatility:
//    here one scale s per expiry row. Along TENOR, at a fixed expiry, the vols are set by how rates
//    of different maturities move together: the mean reversions, the volatility ratio and the
//    correlation. Everything is fitted together by least squares (Levenberg–Marquardt, residuals in
//    bp of normal vol). The one-factor model fitted to the same grid with the same freedom in the
//    level (`factors = 1`: one mean reversion) is the benchmark: the difference in fit is what the
//    second factor buys, and it is all in the tenor direction.
//
//    The fit is not well identified in every direction, which is a property of the model, not of the
//    optimiser: ρ close to −1 with a larger σ_2 gives almost the same ATM matrix as a moderate ρ with
//    a smaller one. Bounds keep the parameters where they mean something; the report says which
//    bounds are active.
//
// 2. LEVEL, from the co-terminal strip of the trade being priced. Holding the shape, the common
//    volatility scale s(t) is bootstrapped one expiry at a time so that every co-terminal swaption
//    is repriced exactly, as ζ is in the one-factor model (lgm_calibration.hpp). The result is a
//    two-factor model that agrees with the one-factor model on every European the Bermudan
//    contains, and differs from it only in how rates on different exercise dates are correlated.

#include "rl/core/date.hpp"
#include "rl/market/curve_set.hpp"
#include "rl/market/swaption_cube.hpp"
#include "rl/math/lm.hpp"
#include "rl/models/g2pp.hpp"
#include "rl/models/lgm_calibration.hpp"

#include <string>
#include <vector>

namespace rl::models {

/// ATM swaptions on the grid expiries × tenors (whole years): the right to enter, in `expiry` years,
/// a SOFR swap of `tenor` years at its forward rate, with the market vol read from the cube.
std::vector<CalibrationInstrument> makeAtmMatrix(core::Date asOf, const market::CurveSet<double>& curves, const market::SwaptionCube<double>& cube,
                                                 const std::vector<int>& expiryYears, const std::vector<int>& tenorYears);

struct G2ppFitOptions {
    int factors = 2;                ///< 1: no second factor, fit κ_1 and the level only
    bool timeDependentScale = true; ///< one volatility scale per expiry row; false: constant parameters
    double kappa1 = 0.30;           ///< starting point
    double kappa2 = 0.02;
    double sigma1 = 0.010;
    double sigma2 = 0.008;
    double rho = -0.70;
    double kappaMin = 1e-3, kappaMax = 1.5;
    double rhoBound = 0.99;         ///< |ρ| ≤ rhoBound
    int maxIterations = 100;
};

struct G2ppFitRow {
    std::string label;
    double expiry = 0.0, tenor = 0.0;
    double marketVolBp = 0.0, modelVolBp = 0.0, errorBp = 0.0;
};

struct G2ppFit {
    G2pp model;
    std::vector<G2ppFitRow> rows;
    double rmsErrorBp = 0.0;
    double maxAbsErrorBp = 0.0;
    int iterations = 0;
    math::LmStatus status = math::LmStatus::Failed;
    std::vector<std::string> activeBounds; ///< parameters that ended on a bound
};

/// Least-squares fit to `instruments`. The shape parameters are κ_1 (and κ_2, σ_2/σ_1, ρ); the level
/// is one scale per distinct expiry, or a single constant. Throws InputError unless there are more
/// instruments than parameters.
G2ppFit fitG2pp(const std::vector<CalibrationInstrument>& instruments, const G2ppFitOptions& options = {});

struct G2ppStripRow {
    std::string label;
    double expiry = 0.0, tenor = 0.0;
    double marketVolBp = 0.0, modelVolBp = 0.0, errorBp = 0.0;
    double scale = 0.0;   ///< s on the interval ending at this expiry
    bool clamped = false; ///< the market price needed less variance than the earlier expiries already imply
};

struct G2ppStripCalibration {
    G2pp model;
    std::vector<G2ppStripRow> rows;
    double maxAbsErrorBp = 0.0;
    bool ok = false; ///< every instrument repriced and nothing was clamped
};

/// Bootstrap the volatility scale of `shape` (its own scale is ignored) to a strip in increasing
/// expiry order. Throws InputError on an empty or unordered strip.
G2ppStripCalibration calibrateG2ppToStrip(const G2pp& shape, const std::vector<CalibrationInstrument>& strip);

/// A constant-parameter model that reprices TWO swaptions exactly, for given mean reversions and
/// correlation: the overall volatility and the ratio σ_2/σ_1 are solved for. With the two swaptions
/// on the two rates of a spread product, models built this way for different ρ agree on both rates'
/// volatilities and differ only in their correlation.
///
/// Not every (κ_1, κ_2, ρ) can do it: the ratio moves the two vols relative to each other only so
/// far, so the pair of vols itself restricts ρ. Throws CalibrationError when no ratio in
/// [ratioMin, ratioMax] works. When more than one does, the one closest to `preferredRatio` is returned.
G2pp calibrateG2ppToPair(double kappa1, double kappa2, double rho, const CalibrationInstrument& first, const CalibrationInstrument& second,
                         double preferredRatio = 1.0, double ratioMin = 0.05, double ratioMax = 20.0);

/// The one-factor counterpart: mean reversion and volatility solved for the same two swaptions.
/// Throws CalibrationError if no mean reversion in [kappaMin, kappaMax] gives the pair.
G2pp calibrateOneFactorToPair(const CalibrationInstrument& first, const CalibrationInstrument& second, double kappaMin = -0.05, double kappaMax = 1.0);

} // namespace rl::models
