// Calibration to the co-terminal strip. The checks are the ones a desk makes every morning: does the
// model reprice what it was calibrated to, and is the fitted variance a variance.
#include "rl/models/lgm_calibration.hpp"

#include "rl/core/error.hpp"
#include "rl/models/lgm_swaption.hpp"
#include "test_support.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>

using namespace rl::models;
using namespace rl::models::test;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

TEST_CASE("Strip from the cube reprices to 1e-6 bp for every mean reversion", "[models][calibration]") {
    const UsdMarket& m = usdMarket();
    for (double strike : {0.025, 0.035, 0.045}) {
        const auto strip = makeStrip(rights10nc1(strike), m.cube);
        REQUIRE(strip.size() == 9);
        for (double kappa : {-0.01, 0.0, 0.03, 0.10}) {
            const LgmCalibration cal = calibrateLgm(kappa, strip);
            INFO("strike " << strike << " kappa " << kappa);
            CHECK(cal.ok);
            CHECK(cal.maxAbsErrorBp < 1e-6);
            double prev = 0.0;
            for (const CalibrationRow& row : cal.rows) {
                CHECK_FALSE(row.clamped);
                CHECK(row.zeta > prev);
                CHECK(row.hullWhiteSigmaBp > 10.0);
                CHECK(row.hullWhiteSigmaBp < 300.0);
                prev = row.zeta;
            }
        }
    }
}

TEST_CASE("Mean reversion is not identified by the strip: prices agree, variances do not", "[models][calibration]") {
    const UsdMarket& m = usdMarket();
    const auto strip = makeStrip(rights10nc1(0.035), m.cube);
    const LgmCalibration lo = calibrateLgm(0.0, strip), hi = calibrateLgm(0.10, strip);
    for (std::size_t i = 0; i < strip.size(); ++i) {
        CHECK_THAT(lgmSwaptionPrice(lo.model, strip[i].right), WithinRel(lgmSwaptionPrice(hi.model, strip[i].right), 1e-9));
        // With mean reversion later cashflows react less to the state (H is concave), so the state
        // needs more variance to produce the same swaption price.
        CHECK(hi.rows[i].zeta > lo.rows[i].zeta);
    }
}

TEST_CASE("Flat normal vol, no mean reversion: the short-rate vol is close to the swaption vol", "[models][calibration]") {
    // With kappa = 0 every rate moves in parallel with the short rate, so a swap rate's normal vol
    // is roughly the short-rate vol: calibrated sigma should sit near the quoted 85 bp.
    const LgmCalibration cal = calibrateLgm(0.0, makeStrip(rights10nc1(0.035), 0.0085));
    CHECK(cal.ok);
    for (const CalibrationRow& row : cal.rows) CHECK_THAT(row.hullWhiteSigmaBp, WithinRel(85.0, 0.10));
}

TEST_CASE("An inconsistent strip is reported, not fitted", "[models][calibration]") {
    // Vol collapsing with expiry so fast that total variance would have to fall.
    auto strip = makeStrip(rights10nc1(0.035), 0.0085);
    for (std::size_t i = 4; i < strip.size(); ++i) {
        strip[i].marketVol = 0.0010;
        strip[i].marketPrice = marketSwaptionPrice(strip[i].right, strip[i].marketVol);
    }
    const LgmCalibration cal = calibrateLgm(0.03, strip);
    CHECK_FALSE(cal.ok);
    CHECK(cal.rows[4].clamped);
    CHECK(cal.rows[4].errorBp > 1.0); // the model is left visibly too high there
    CHECK(cal.rows[3].clamped == false);
    CHECK(std::fabs(cal.rows[3].errorBp) < 1e-6);
}

TEST_CASE("Calibration input validation", "[models][calibration]") {
    CHECK_THROWS_AS(calibrateLgm(0.03, {}), rl::InputError);
    auto strip = makeStrip(rights10nc1(0.035), 0.0085);
    std::swap(strip[0], strip[1]);
    CHECK_THROWS_AS(calibrateLgm(0.03, strip), rl::InputError);
}
