#include "rl/models/lgm_calibration.hpp"

#include "rl/core/error.hpp"
#include "rl/math/rootfind.hpp"
#include "rl/models/lgm_swaption.hpp"
#include "rl/vol/cube.hpp"

#include <algorithm>
#include <cmath>

namespace rl::models {

std::vector<CalibrationInstrument> makeStrip(const std::vector<ExerciseRight>& rights, const market::SwaptionCube<double>& cube) {
    std::vector<CalibrationInstrument> strip;
    for (const ExerciseRight& r : rights) {
        CalibrationInstrument ci;
        ci.right = r;
        ci.marketVol = vol::cubeNormalVol<double>(cube, r.time, r.tenorYears, r.strike, r.forwardRate());
        ci.marketPrice = marketSwaptionPrice(r, ci.marketVol);
        strip.push_back(std::move(ci));
    }
    return strip;
}

std::vector<CalibrationInstrument> makeStrip(const std::vector<ExerciseRight>& rights, double flatNormalVol) {
    std::vector<CalibrationInstrument> strip;
    for (const ExerciseRight& r : rights) strip.push_back({r, flatNormalVol, marketSwaptionPrice(r, flatNormalVol)});
    return strip;
}

LgmCalibration calibrateLgm(double kappa, const std::vector<CalibrationInstrument>& strip) {
    RL_REQUIRE(!strip.empty(), "calibrateLgm: empty strip");
    for (std::size_t i = 1; i < strip.size(); ++i)
        RL_REQUIRE(strip[i].right.time > strip[i - 1].right.time, "calibrateLgm: the strip must be in increasing expiry order");

    // Only H is needed while solving; any valid ζ node will do for the helper model.
    const Lgm1F shape(kappa, {1.0}, {1.0});
    std::vector<double> times, zetas;
    std::vector<bool> clamped;
    double previous = 0.0;
    for (const CalibrationInstrument& ci : strip) {
        auto error = [&](double logZeta) { return lgmOptionValue(shape, std::exp(logZeta), ci.right.flows) - ci.marketPrice; };
        // Solve in log ζ: ζ spans orders of magnitude (≈ σ² t, 1e-6 … 1e-2) and must stay positive.
        double lo = std::log(1e-12), hi = std::log(1.0);
        double zeta;
        if (error(lo) >= 0.0) {
            zeta = std::exp(lo); // the market price is at or below intrinsic: no variance needed
        } else {
            RL_REQUIRE(error(hi) > 0.0, "calibrateLgm: market price of " + ci.right.label + " is not attainable");
            math::RootOptions ro;
            ro.xTol = 1e-13;
            ro.maxIter = 200;
            zeta = std::exp(math::brent(error, lo, hi, ro).x);
        }
        const bool clamp = zeta < previous;
        if (clamp) zeta = previous;
        // Keep the nodes strictly usable by Lgm1F even when clamped.
        times.push_back(ci.right.time);
        zetas.push_back(zeta);
        clamped.push_back(clamp);
        previous = zeta;
    }

    LgmCalibration out{Lgm1F(kappa, times, zetas), {}, 0.0, true};
    const std::vector<double> sigmas = out.model.hullWhiteSigmas();
    for (std::size_t i = 0; i < strip.size(); ++i) {
        const ExerciseRight& r = strip[i].right;
        CalibrationRow row;
        row.label = r.label;
        row.expiry = r.time;
        row.tenor = r.tenorYears;
        row.forward = r.forwardRate();
        row.strike = r.strike;
        row.marketVolBp = strip[i].marketVol * 1e4;
        row.modelVolBp = impliedNormalVol(r, lgmSwaptionPrice(out.model, r)) * 1e4;
        row.errorBp = row.modelVolBp - row.marketVolBp;
        row.zeta = zetas[i];
        row.hullWhiteSigmaBp = sigmas[i] * 1e4;
        row.clamped = clamped[i];
        out.maxAbsErrorBp = std::max(out.maxAbsErrorBp, std::fabs(row.errorBp));
        if (row.clamped) out.ok = false;
        out.rows.push_back(std::move(row));
    }
    if (out.maxAbsErrorBp > 0.01) out.ok = false;
    return out;
}

} // namespace rl::models
