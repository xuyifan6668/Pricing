// A swap's expected exposure is a strip of European swaptions: closed form against simulation.
#include "rl/pricers/exposure.hpp"

#include "rl/linear/leg_pv.hpp"
#include "rl/models/lgm_calibration.hpp"
#include "test_support.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>

using namespace rl::pricers;
using namespace rl::pricers::test;
using Catch::Matchers::WithinAbs;

TEST_CASE("Expected exposure: closed-form swaption strip matches Monte Carlo", "[pricers][exposure]") {
    const auto& m = usdMarket();
    const auto swap = spotSwap(10, 0.035);
    const CalibratedBermudan b = calibrateBermudan(m, swap, 1.0, 0, 0.03);
    const auto profile = swapExposureProfile(b.problem.model(), swap, m.curves);
    REQUIRE(profile.size() == 10); // the spot start, then the nine later period starts
    const auto european = b.problem.europeanValues();
    for (std::size_t i = 0; i < profile.size(); ++i) {
        const ExposurePoint& p = profile[i];
        CHECK_THAT(p.expectedExposureMc, WithinAbs(p.expectedExposure, 4.0 * p.mcStandardError));
        // From year 1 on these are exactly the Bermudan's calibration swaptions, so EE is also the MARKET swaption price.
        if (i > 0) CHECK_THAT(p.expectedExposure, WithinAbs(european[i - 1], 1e-14));
        CHECK(p.potentialExposure > p.expectedExposure);
        CHECK(p.expectedNegativeExposure <= 0.0);
    }
}

TEST_CASE("Payer and receiver exposures: EE(payer) + ENE(receiver) = 0", "[pricers][exposure]") {
    const auto& m = usdMarket();
    const auto pay = spotSwap(10, 0.035, true), rec = spotSwap(10, 0.035, false);
    const rl::models::Lgm1F model = rl::models::Lgm1F::constantSigma(0.03, 0.009, 30.0);
    const auto pp = swapExposureProfile(model, pay, m.curves), pr = swapExposureProfile(model, rec, m.curves);
    for (std::size_t i = 0; i < pp.size(); ++i) {
        CHECK_THAT(pp[i].expectedExposure + pr[i].expectedNegativeExposure, WithinAbs(0.0, 1e-14));
        CHECK_THAT(pr[i].expectedExposure + pp[i].expectedNegativeExposure, WithinAbs(0.0, 1e-14));
    }
}

TEST_CASE("Exposure of a par swap rises then falls: diffusion against amortisation", "[pricers][exposure]") {
    // Variance of rates grows with time, the remaining annuity shrinks: the profile is hump-shaped.
    const auto& m = usdMarket();
    const auto swap = spotSwap(10, 0.0345);
    const auto profile = swapExposureProfile(rl::models::Lgm1F::constantSigma(0.03, 0.009, 30.0), swap, m.curves);
    std::size_t peak = 0;
    for (std::size_t i = 0; i < profile.size(); ++i)
        if (profile[i].expectedExposure > profile[peak].expectedExposure) peak = i;
    CHECK(peak > 0);
    CHECK(peak + 1 < profile.size());
    CHECK(profile.back().expectedExposure < profile[peak].expectedExposure);
}
