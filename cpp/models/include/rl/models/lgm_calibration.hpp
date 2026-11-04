#pragma once
// Calibrating LGM to a co-terminal swaption strip.
//
// Which swaptions? A Bermudan that can be exercised at t_1 < … < t_n into the swap ending at T is,
// at each date, a choice between exercising a European into the remaining swap and waiting. The
// model must at least price each of those Europeans correctly: the i-th has expiry t_i and
// underlying [t_i, T]. They all end on the same date, hence "co-terminal". Calibrating instead to
// swaptions of one fixed tenor ("co-initial") would price a different set of options and leave the
// Bermudan's own exercise values wrong.
//
// At which strike? The Bermudan's own fixed rate, read off the smile of the vol cube. A one-factor
// Gaussian model has no smile of its own to speak of, so it is calibrated where the trade lives.
//
// How? With H (that is, κ) fixed, the i-th swaption's model price depends on the model only through
// ζ(t_i), and increases with it. So ζ(t_1), ζ(t_2), … are found one at a time by 1-D root search.
// The result must be non-decreasing in time, since ζ is accumulated variance. If the market strip
// asks for a DECREASE, the vols are inconsistent with a one-factor model at this κ; that is
// reported per instrument and ζ is held flat there, leaving a visible repricing error rather than
// a silently "fitted" model.

#include "rl/market/swaption_cube.hpp"
#include "rl/models/flows.hpp"
#include "rl/models/lgm.hpp"

#include <string>
#include <vector>

namespace rl::models {

struct CalibrationInstrument {
    ExerciseRight right;
    double marketVol = 0.0;   ///< normal vol, decimal
    double marketPrice = 0.0; ///< annuity × Bachelier
};

struct CalibrationRow {
    std::string label;
    double expiry = 0.0;        ///< years
    double tenor = 0.0;         ///< years
    double forward = 0.0;
    double strike = 0.0;
    double marketVolBp = 0.0;
    double modelVolBp = 0.0;
    double errorBp = 0.0;       ///< model − market, bp of normal vol
    double zeta = 0.0;
    double hullWhiteSigmaBp = 0.0; ///< short-rate vol on the interval ending at this expiry, bp
    bool clamped = false;       ///< the market asked for decreasing ζ here
};

struct LgmCalibration {
    Lgm1F model;
    std::vector<CalibrationRow> rows;
    double maxAbsErrorBp = 0.0;
    bool ok = false; ///< every instrument repriced and nothing was clamped
};

/// The strip for a set of exercise rights: market vol from the cube at each right's expiry, tenor,
/// strike and forward.
std::vector<CalibrationInstrument> makeStrip(const std::vector<ExerciseRight>& rights, const market::SwaptionCube<double>& cube);

/// The strip from flat normal vol (testing, and "what if the surface were flat").
std::vector<CalibrationInstrument> makeStrip(const std::vector<ExerciseRight>& rights, double flatNormalVol);

/// Bootstrap ζ at the strip's expiries for the given mean reversion. Instruments must be in
/// increasing expiry order. Throws InputError on an empty or unordered strip.
LgmCalibration calibrateLgm(double kappa, const std::vector<CalibrationInstrument>& strip);

/// Marking the mean reversion.
///
/// A co-terminal strip pins ζ for ANY κ, so it cannot say what κ is. A second set of swaptions with
/// the same expiries but a different tenor can: given ζ(t), the model vol of a short swap starting
/// at t scales with H(t + tenor) − H(t) = e^{−κt} H(tenor), which falls with κ, while the co-terminal
/// swap's does not move by construction. So κ is chosen to minimise the RMS error on the second set,
/// re-bootstrapping ζ to the strip at every trial value. The second set is fitted in the
/// least-squares sense only: one number cannot reprice nine instruments, and the residual says how
/// far a one-factor model is from holding both strips at once.
struct MeanReversionCalibration {
    double kappa = 0.0;
    LgmCalibration calibration;        ///< the strip calibration at that κ: still exact
    double secondaryRmsBp = 0.0;       ///< RMS error on the second set, bp of normal vol
    std::vector<double> secondaryErrorsBp; ///< model − market per instrument of the second set
    std::vector<double> scanKappa, scanRmsBp; ///< the objective over the search range
    bool atBound = false;              ///< the minimum is at an end of the range: κ is not identified inside it
};

/// Throws InputError on an empty second set or an empty range.
MeanReversionCalibration calibrateLgmMeanReversion(const std::vector<CalibrationInstrument>& strip, const std::vector<CalibrationInstrument>& secondary,
                                                   double kappaMin = -0.02, double kappaMax = 0.25);

} // namespace rl::models
