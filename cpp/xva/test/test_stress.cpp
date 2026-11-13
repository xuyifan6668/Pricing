// Stress: the whole chain re-run in shocked markets, checked for the same properties as in the base.
//
// Every other test runs on one calm market. A stress suite asks whether the chain still holds
// together far from it: the curve re-bootstraps, the model re-calibrates exactly, the Bermudan stays
// inside its Monte Carlo bracket and above its Europeans, and the exposure engine still reproduces
// the closed form. And where the chain cannot work, it must fail cleanly, not return a number.
#include "rl/core/error.hpp"
#include "rl/pricers/bermudan_mc.hpp"
#include "rl/pricers/bermudan_pde.hpp"
#include "rl/pricers/hedging.hpp"
#include "rl/xva/exposure_engine.hpp"
#include "rl/xva/xva.hpp"
#include "test_support.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <algorithm>
#include <cmath>

using namespace rl;
using namespace rl::xva;
using namespace rl::xva::test;
using Catch::Matchers::WithinAbs;

namespace {
const double kBp = 1e-4;

pricers::MarketMove move(double parallel, double slope, double vol) { return pricers::ratesMove(usdMarket(), parallel, slope, vol); }
} // namespace

TEST_CASE("Stress: the chain holds together in shocked markets", "[xva][stress]") {
    const struct {
        const char* name;
        double parallel, slope, vol;
    } scenarios[] = {{"base", 0, 0, 0},
                     {"rates +200 bp", 200, 0, 0},
                     {"rates -250 bp", -250, 0, 0},
                     {"inverted: short end +150, long end -100", 25, -250, 0},
                     {"vol +60 bp", 0, 0, 60},
                     {"vol -35 bp", 0, 0, -35},
                     {"rates +150 bp and vol +40 bp", 150, 0, 40},
                     {"rates +500 bp", 500, 0, 0}};
    for (const auto& sc : scenarios) {
        INFO(sc.name);
        const models::UsdMarket stressed = pricers::applyMove(usdMarket(), move(sc.parallel, sc.slope, sc.vol));
        const auto& disc = stressed.curves.discount();
        CHECK(disc.df(1.0) > 0.0);
        CHECK(disc.df(30.0) < disc.df(1.0));
        // A swap struck at the stressed market's own par rate.
        const instruments::Swap probe = spotSwap(10, 0.035);
        const double strike = models::makeExerciseRight(probe, stressed.curves, std::get<instruments::FixedCoupon>(probe.fixedLeg().cashflows().front()).accrualStart).forwardRate();
        const instruments::Swap trade = spotSwap(10, strike);
        const pricers::CalibratedBermudan b = pricers::calibrateBermudan(stressed, trade, 1.0, 0, 0.03);
        CHECK(b.calibration.ok);
        CHECK(b.calibration.maxAbsErrorBp < 1e-6);

        const double pde = pricers::priceBermudanPde(b.problem).value;
        const std::vector<double> eu = b.problem.europeanValues();
        double total = 0.0;
        for (double e : eu) total += e;
        CHECK(pde >= *std::max_element(eu.begin(), eu.end()));
        CHECK(pde <= total);
        const pricers::LsmcResult lo = pricers::priceBermudanLsmc(b.problem, {.trainingPaths = 10000, .pricingPaths = 40000});
        const pricers::UpperBoundResult up = pricers::priceBermudanUpperBound(b.problem, lo, {.outerPaths = 300, .innerPaths = 300});
        INFO("lower " << lo.lower / kBp << ", PDE " << pde / kBp << ", upper " << up.upper / kBp << " bp");
        CHECK(lo.lower <= pde + 4.0 * lo.standardError + 0.1 * kBp);
        CHECK(up.upper >= pde - 4.0 * up.standardError - 0.1 * kBp);

        // The exposure engine under the stressed model still reproduces the swaption strip.
        const SwapExposure exact = swapExposureClosedForm(b.problem.model(), trade, stressed.curves);
        const ExposureProfile mc = simulateExposure(b.problem.model(), {trade}, stressed.curves, exact.times, std::nullopt, {.paths = 20000});
        for (std::size_t i = 1; i < exact.times.size(); ++i) CHECK_THAT(mc.epe[i], WithinAbs(exact.epe[i], 4.0 * mc.epeStandardError[i]));
    }
}

TEST_CASE("Stress: a scenario that moves basis spreads like rates is not solvable, and the bootstrap says so", "[xva][stress]") {
    // Shifting EVERY quote by +200 bp also moves the Fed Funds-SOFR basis from -4 bp to +196 bp. The
    // first version of these scenarios did that by mistake; the bootstrap refused rather than return a curve.
    pricers::MarketMove all;
    all.quoteShiftBp.assign(usdMarket().problem.size(), 200.0);
    CHECK_THROWS_AS(pricers::applyMove(usdMarket(), all), rl::Error);
}

TEST_CASE("Stress: below zero rates the unshifted cube cannot be read, and the chain says so", "[xva][stress]") {
    // The cube's SABR has beta 0.5 and no shift, so it is undefined for a non-positive forward. With
    // rates 450 bp lower the 10-year forward is negative: calibration must refuse, not return NaN.
    const models::UsdMarket stressed = pricers::applyMove(usdMarket(), move(-450, 0, 0));
    const instruments::Swap trade = spotSwap(10, 0.0);
    CHECK_THROWS_AS(pricers::calibrateBermudan(stressed, trade, 1.0, 0, 0.03), rl::Error);
}
