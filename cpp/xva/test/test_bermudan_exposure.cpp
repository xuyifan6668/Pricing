// The Bermudan's exposure: the same simulation must reproduce its price, and its limits are known.
#include "rl/xva/bermudan_exposure.hpp"

#include "rl/core/error.hpp"
#include "rl/models/g2pp_calibration.hpp"
#include "rl/pricers/bermudan_g2pp.hpp"
#include "rl/models/lgm_swaption.hpp"
#include "rl/pricers/bermudan_pde.hpp"
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
pricers::CalibratedBermudan bermudan(double rate, bool payer = true) { return pricers::calibrateBermudan(usdMarket(), spotSwap(10, rate, payer), 1.0, 0, 0.03); }
const double kBp = 1e-4;
} // namespace

TEST_CASE("Bermudan exposure: the grid and the simulation both reproduce the PDE price", "[xva][bermudan]") {
    for (bool payer : {true, false}) {
        const pricers::CalibratedBermudan b = bermudan(0.035, payer);
        const double pde = pricers::priceBermudanPde(b.problem, {.spaceSteps = 800, .timeSteps = 800}).value;
        const BermudanExposure e = bermudanExposure(b.problem);
        CHECK_THAT(e.price, WithinAbs(pde, 0.05 * kBp));
        CHECK_THAT(e.simulatedPrice, WithinAbs(e.price, 4.0 * e.simulatedPriceStandardError));
        CHECK(e.simulatedPriceStandardError < 3.0 * kBp);
    }
}

TEST_CASE("Bermudan exposure: an option until exercised", "[xva][bermudan]") {
    const pricers::CalibratedBermudan b = bermudan(0.035);
    const BermudanExposure e = bermudanExposure(b.problem);
    REQUIRE(e.times.size() == 10);
    CHECK(e.times[0] == 0.0);
    CHECK_THAT(e.epe[0], WithinAbs(e.price, 1e-15));
    // On the first exercise date the position is worth max(exercise, continuation) on every path: a
    // martingale from today, so its mean is the price, and it cannot be negative.
    CHECK_THAT(e.epe[1], WithinAbs(e.price, 4.0 * e.epeStandardError[1]));
    CHECK(e.ene[1] == 0.0);
    for (std::size_t i = 1; i < e.times.size(); ++i) {
        CHECK(e.epe[i] >= e.optionEpe[i]);
        CHECK(e.ene[i] <= 0.0);
        CHECK(e.exercisedBy[i] >= e.exercisedBy[i - 1]); // exercise is once and for all
        CHECK(e.exercisedBy[i] <= 1.0);
    }
    CHECK(e.exercisedBy.back() > 0.4);
    CHECK(e.optionEpe.back() == 0.0);  // nothing left to hold after the last date
    CHECK(e.ene.back() < -1.0 * kBp);  // swaps entered earlier have since gone against the holder
    CHECK(e.epe.back() < 0.5 * e.epe[1]);
}

TEST_CASE("Bermudan exposure: one exercise date is the European; deep in the money is the swap", "[xva][bermudan]") {
    const pricers::CalibratedBermudan b = bermudan(0.035);
    const pricers::BermudanProblem single(b.problem.model(), {b.problem.rights()[3]});
    const BermudanExposure e = bermudanExposure(single);
    CHECK_THAT(e.price, WithinAbs(models::lgmSwaptionPrice(single.model(), single.rights()[0]), 0.01 * kBp));
    CHECK_THAT(e.epe[1], WithinAbs(e.price, 4.0 * e.epeStandardError[1]));

    // A payer struck far below the market is exercised at once on almost every path. From then on the holder
    // owns the swap, whose exposure is the closed-form swaption strip.
    const double strike = 0.001;
    const pricers::BermudanProblem deep(b.problem.model(), pricers::calibrateBermudan(usdMarket(), spotSwap(10, strike), 1.0, 0, 0.03).problem.rights());
    const BermudanExposure d = bermudanExposure(deep);
    const SwapExposure swap = swapExposureClosedForm(deep.model(), spotSwap(10, strike), usdMarket().curves);
    CHECK(d.exercisedBy[1] > 0.98);
    for (std::size_t i = 1; i < d.times.size(); ++i) {
        // swap.times is 0, the spot date, then the anniversaries: date i of the Bermudan is entry i + 1.
        CHECK_THAT(d.times[i], WithinAbs(swap.times[i + 1], 1e-9));
        CHECK_THAT(d.epe[i], WithinAbs(swap.epe[i + 1], 4.0 * d.epeStandardError[i] + 1e-6));
    }
}

TEST_CASE("Bermudan exposure: more rights to exercise mean more exposure at the start", "[xva][bermudan]") {
    const pricers::CalibratedBermudan b = bermudan(0.035);
    const BermudanExposure all = bermudanExposure(b.problem);
    const pricers::BermudanProblem late(b.problem.model(), {b.problem.rights()[6], b.problem.rights()[7], b.problem.rights()[8]});
    CHECK(all.price > bermudanExposure(late).price);
    CHECK_THROWS_AS(bermudanExposure(b.problem, {.paths = 10}), rl::InputError);
    CHECK_THROWS_AS(bermudanExposure(b.problem, {.gridPoints = 100}), rl::InputError);
}

TEST_CASE("Two-factor Bermudan exposure with the second factor off matches the one-factor grid", "[xva][bermudan][g2pp]") {
    const pricers::CalibratedBermudan one = bermudan(0.035);
    const BermudanExposure grid = bermudanExposure(one.problem, {.paths = 200000});
    const pricers::CalibratedBermudanG2pp two = pricers::calibrateBermudanG2pp(usdMarket(), spotSwap(10, 0.035), 1.0, 0, models::G2pp(0.03, 1.0, 0.01, 0.0, 0.0));
    const BermudanExposure amc = bermudanExposureG2pp(two.calibration.model, two.rights);
    // A regression rule is a lower bound, and a good one: within a basis point of the optimal price.
    CHECK(amc.price <= grid.price + 4.0 * amc.simulatedPriceStandardError);
    CHECK(amc.price >= grid.price - 4.0 * amc.simulatedPriceStandardError - 1.0 * kBp);
    REQUIRE(amc.times.size() == grid.times.size());
    for (std::size_t i = 1; i < grid.times.size(); ++i) {
        INFO("date " << i << ": regression " << amc.epe[i] / kBp << " bp, grid " << grid.epe[i] / kBp << " bp");
        CHECK_THAT(amc.epe[i], WithinAbs(grid.epe[i], 0.01 * grid.epe[i] + 4.0 * (amc.epeStandardError[i] + grid.epeStandardError[i])));
        // The split between option and swap is a price too, so it must agree. (Exercise frequencies need
        // not: the two simulations run under different numeraires.)
        CHECK_THAT(amc.optionEpe[i], WithinAbs(grid.optionEpe[i], 0.02 * grid.epe[i] + 4.0 * (amc.epeStandardError[i] + grid.epeStandardError[i])));
    }
}

TEST_CASE("Two-factor Bermudan exposure: an option until exercised, whatever the correlation", "[xva][bermudan][g2pp]") {
    for (double rho : {-0.7, 0.5}) {
        const pricers::CalibratedBermudanG2pp b = pricers::calibrateBermudanG2pp(usdMarket(), spotSwap(10, 0.035), 1.0, 0, models::G2pp(0.40, 0.03, 0.012, 0.009, rho));
        const BermudanExposure e = bermudanExposureG2pp(b.calibration.model, b.rights);
        // The independent pricing regression gives the same lower bound, to sampling error.
        const pricers::G2ppLsmcResult ref = pricers::priceBermudanLsmcG2pp(b.calibration.model, b.rights);
        CHECK_THAT(e.price, WithinAbs(ref.lower, 4.0 * (e.simulatedPriceStandardError + ref.standardError) + 1.0 * kBp));
        CHECK_THAT(e.epe[1], WithinAbs(e.price, 0.01 * e.price + 4.0 * e.epeStandardError[1])); // a martingale to the first date
        CHECK(e.ene[1] == 0.0);
        for (std::size_t i = 1; i < e.times.size(); ++i) {
            CHECK(e.epe[i] >= e.optionEpe[i]);
            CHECK(e.exercisedBy[i] >= e.exercisedBy[i - 1]);
        }
        CHECK(e.optionEpe.back() == 0.0);
    }
    const pricers::CalibratedBermudanG2pp b = pricers::calibrateBermudanG2pp(usdMarket(), spotSwap(10, 0.035), 1.0, 0, models::G2pp(0.40, 0.03, 0.012, 0.009, -0.7));
    CHECK_THROWS_AS(bermudanExposureG2pp(b.calibration.model, {}), rl::InputError);
    CHECK_THROWS_AS(bermudanExposureG2pp(b.calibration.model, {b.rights[2], b.rights[1]}), rl::InputError);
}
