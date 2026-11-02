// G2++ calibration: recover known parameters, beat the one-factor fit, reprice a strip exactly.
#include "rl/models/g2pp_calibration.hpp"

#include "rl/core/error.hpp"
#include "rl/models/g2pp_swaption.hpp"
#include "rl/models/lgm_swaption.hpp"
#include "test_support.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>

using namespace rl::models;
using namespace rl::models::test;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace {
std::vector<CalibrationInstrument> matrix() {
    const UsdMarket& m = usdMarket();
    return makeAtmMatrix(m.asOf, m.curves, m.cube, {1, 2, 5, 10}, {1, 2, 5, 10});
}
/// The same instruments with market vols replaced by a model's own.
std::vector<CalibrationInstrument> repricedBy(const G2pp& model, std::vector<CalibrationInstrument> ins) {
    for (CalibrationInstrument& ci : ins) {
        ci.marketPrice = g2ppSwaptionPrice(model, ci.right);
        ci.marketVol = impliedNormalVol(ci.right, ci.marketPrice);
    }
    return ins;
}
} // namespace

TEST_CASE("ATM matrix instruments are at the money, with the cube's ATM vol", "[models][g2pp][calibration]") {
    const auto ins = matrix();
    REQUIRE(ins.size() == 16);
    for (const CalibrationInstrument& ci : ins) {
        CHECK_THAT(ci.right.forwardValue(), WithinAbs(0.0, 1e-12));
        CHECK_THAT(ci.right.forwardRate(), WithinAbs(ci.right.strike, 1e-12));
        CHECK(ci.marketVol > 0.004);
        CHECK(ci.marketVol < 0.015);
    }
    CHECK(ins.front().right.label == "1Yx1Y");
    CHECK_THAT(ins.back().right.time, WithinAbs(10.0, 0.05));
    CHECK_THAT(ins.back().right.tenorYears, WithinAbs(10.0, 0.05));
}

TEST_CASE("The fit recovers a two-factor model's own vols", "[models][g2pp][calibration]") {
    const G2pp truth(0.40, 0.03, 0.012, 0.009, -0.6);
    const G2ppFit fit = fitG2pp(repricedBy(truth, matrix()));
    CHECK(fit.maxAbsErrorBp < 1e-2);
    // The vols are what must match; the correlation they imply comes back with them.
    CHECK_THAT(fit.model.forwardRateCorrelation(1.0, 10.0), WithinAbs(truth.forwardRateCorrelation(1.0, 10.0), 0.05));
    REQUIRE(fit.model.scales().size() == 4); // one per expiry row
    for (double s : fit.model.scales()) CHECK_THAT(s * fit.model.sigma1(), WithinRel(truth.sigma1(), 0.05));
}

TEST_CASE("One factor fitted to a two-factor matrix leaves an error; the one-factor fit is exact on one-factor data", "[models][g2pp][calibration]") {
    const G2pp truth(0.40, 0.03, 0.012, 0.009, -0.6);
    const auto ins = repricedBy(truth, matrix());
    G2ppFitOptions one;
    one.factors = 1;
    one.kappa1 = 0.05;
    const G2ppFit fit1 = fitG2pp(ins, one);
    CHECK(fit1.model.sigma2() == 0.0);
    CHECK(fit1.rmsErrorBp > 1.0); // a single mean reversion cannot bend the matrix in both directions

    const G2pp flat(0.06, 1.0, 0.009, 0.0, 0.0);
    const G2ppFit exact = fitG2pp(repricedBy(flat, matrix()), one);
    CHECK(exact.maxAbsErrorBp < 1e-4);
    CHECK_THAT(exact.model.kappa1(), WithinAbs(0.06, 1e-5));
    for (double s : exact.model.scales()) CHECK_THAT(s * exact.model.sigma1(), WithinRel(0.009, 1e-5));
    // With constant parameters the same data is recovered by two numbers.
    one.timeDependentScale = false;
    const G2ppFit constant = fitG2pp(repricedBy(flat, matrix()), one);
    CHECK(constant.maxAbsErrorBp < 1e-4);
    CHECK(constant.model.scales().empty());
    CHECK_THAT(constant.model.sigma1(), WithinRel(0.009, 1e-5));
    CHECK_THROWS_AS(fitG2pp({ins[0], ins[1]}, G2ppFitOptions{}), rl::InputError);
}

TEST_CASE("On the market matrix two factors fit at least as well as one", "[models][g2pp][calibration]") {
    const auto ins = matrix();
    G2ppFitOptions one;
    one.factors = 1;
    one.kappa1 = 0.05;
    const G2ppFit fit1 = fitG2pp(ins, one), fit2 = fitG2pp(ins);
    CHECK(fit2.rmsErrorBp <= fit1.rmsErrorBp + 1e-6);
    CHECK(fit2.rows.size() == ins.size());
}

TEST_CASE("The volatility scale reprices a co-terminal strip exactly, whatever the shape", "[models][g2pp][calibration]") {
    const auto strip = makeStrip(rights10nc1(0.035), usdMarket().cube);
    for (const G2pp& shape : {G2pp(0.40, 0.03, 0.012, 0.009, -0.6), G2pp(0.03, 1.0, 0.01, 0.0, 0.0), G2pp(0.9, 0.01, 0.02, 0.012, -0.95)}) {
        const G2ppStripCalibration cal = calibrateG2ppToStrip(shape, strip);
        CHECK(cal.ok);
        CHECK(cal.maxAbsErrorBp < 1e-6);
        REQUIRE(cal.rows.size() == strip.size());
        for (const G2ppStripRow& row : cal.rows) CHECK(row.scale > 0.0);
        // Shape is untouched: the correlation structure is the one chosen.
        CHECK_THAT(cal.model.forwardRateCorrelation(1.0, 10.0), WithinAbs(shape.forwardRateCorrelation(1.0, 10.0), 1e-15));
    }
}

TEST_CASE("With no second factor the strip calibration is the one-factor calibration", "[models][g2pp][calibration]") {
    const auto strip = makeStrip(rights10nc1(0.035), usdMarket().cube);
    const double kappa = 0.03, sigma = 0.01;
    const G2ppStripCalibration two = calibrateG2ppToStrip(G2pp(kappa, 1.0, sigma, 0.0, 0.0), strip);
    const LgmCalibration one = calibrateLgm(kappa, strip);
    const std::vector<double> hw = one.model.hullWhiteSigmas();
    for (std::size_t i = 0; i < strip.size(); ++i) CHECK_THAT(two.rows[i].scale * sigma, WithinRel(hw[i], 1e-7));
}

TEST_CASE("An inconsistent strip is flagged, not silently fitted", "[models][g2pp][calibration]") {
    auto strip = makeStrip(rights10nc1(0.035), 0.0080);
    strip[4].marketVol = 0.0020; // far less variance than the first four expiries already imply
    strip[4].marketPrice = marketSwaptionPrice(strip[4].right, strip[4].marketVol);
    const G2ppStripCalibration cal = calibrateG2ppToStrip(G2pp(0.40, 0.03, 0.012, 0.009, -0.6), strip);
    CHECK_FALSE(cal.ok);
    CHECK(cal.rows[4].clamped);
    CHECK(cal.rows[4].errorBp > 1.0);
    CHECK_THROWS_AS(calibrateG2ppToStrip(G2pp(0.1, 0.1, 0.01, 0.01, 0.0), {}), rl::InputError);
}
