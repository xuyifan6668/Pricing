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

} // namespace rl::models
