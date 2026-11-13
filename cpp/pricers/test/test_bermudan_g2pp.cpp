// The two-factor Bermudan, tested on what must hold whatever the parameters, and against the
// one-factor PDE in the limit where the two models are the same model.
#include "rl/pricers/bermudan_g2pp.hpp"

#include "rl/core/error.hpp"
#include "rl/pricers/bermudan_pde.hpp"
#include "test_support.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <algorithm>
#include <cmath>

using namespace rl;
using namespace rl::pricers;
using namespace rl::pricers::test;
using Catch::Matchers::WithinAbs;

namespace {
const models::G2pp kShape(0.40, 0.03, 0.012, 0.009, -0.6);
const double kBp = 1e-4;
} // namespace

TEST_CASE("G2++ strip calibration reprices the Bermudan's own Europeans", "[pricers][g2pp]") {
    const CalibratedBermudanG2pp b = calibrateBermudanG2pp(usdMarket(), spotSwap(10, 0.035), 1.0, 0, kShape);
    REQUIRE(b.rights.size() == 9);
    CHECK(b.calibration.ok);
    CHECK(b.calibration.maxAbsErrorBp < 1e-6);
    // So the two-factor model and the one-factor model agree on every European in the trade.
    const CalibratedBermudan one = bermudan10nc1();
    const std::vector<double> ref = one.problem.europeanValues();
    const G2ppLsmcResult r = priceBermudanLsmcG2pp(b.calibration.model, b.rights, {.trainingPaths = 2000, .pricingPaths = 2000});
    for (std::size_t i = 0; i < ref.size(); ++i) CHECK_THAT(r.europeans[i], WithinAbs(ref[i], 1e-9));
}

TEST_CASE("G2++ Bermudan with one exercise date is the European", "[pricers][g2pp]") {
    const CalibratedBermudanG2pp b = calibrateBermudanG2pp(usdMarket(), spotSwap(10, 0.035), 1.0, 0, kShape);
    for (std::size_t i : {std::size_t(0), std::size_t(4), std::size_t(8)}) {
        const G2ppLsmcResult r = priceBermudanLsmcG2pp(b.calibration.model, {b.rights[i]});
        CHECK_THAT(r.lower, WithinAbs(r.europeans[0], 4.0 * r.standardError));
        CHECK(r.standardError < 2.0 * kBp);
    }
}

TEST_CASE("G2++ Bermudan with the second factor switched off matches the one-factor PDE", "[pricers][g2pp]") {
    for (bool payer : {true, false}) {
        const CalibratedBermudan one = bermudan10nc1(0.035, payer, 0.03);
        const double pde = priceBermudanPde(one.problem).value;
        const CalibratedBermudanG2pp b = calibrateBermudanG2pp(usdMarket(), spotSwap(10, 0.035, payer), 1.0, 0, models::G2pp(0.03, 1.0, 0.01, 0.0, 0.0));
        const G2ppLsmcResult r = priceBermudanLsmcG2pp(b.calibration.model, b.rights);
        // A lower bound: at most the PDE up to sampling error, and not far below it (the rule is nearly optimal).
        CHECK(r.lower <= pde + 4.0 * r.standardError);
        CHECK(r.lower >= pde - 4.0 * r.standardError - 0.5 * kBp);
        CHECK(r.standardError < 1.5 * kBp);
    }
}

TEST_CASE("G2++ Bermudan: at least the best European, at most their sum; in-sample is biased high", "[pricers][g2pp]") {
    for (double strike : {0.025, 0.035, 0.045}) {
        const CalibratedBermudanG2pp b = calibrateBermudanG2pp(usdMarket(), spotSwap(10, strike), 1.0, 0, kShape);
        const G2ppLsmcResult r = priceBermudanLsmcG2pp(b.calibration.model, b.rights);
        const double best = *std::max_element(r.europeans.begin(), r.europeans.end());
        double total = 0.0, probability = 0.0;
        for (double e : r.europeans) total += e;
        for (double p : r.exerciseProbability) probability += p;
        CHECK(r.lower > best + 10.0 * kBp); // the right to choose the date is worth something
        CHECK(r.lower < total);
        CHECK(r.inSample > r.lower - 4.0 * r.standardError);
        CHECK(probability <= 1.0 + 1e-12);
        CHECK(probability > 0.2);
    }
}

TEST_CASE("G2++ Bermudan: more exercise rights are worth more; the run is reproducible", "[pricers][g2pp]") {
    const CalibratedBermudanG2pp b = calibrateBermudanG2pp(usdMarket(), spotSwap(10, 0.035), 1.0, 0, kShape);
    const G2ppLsmcResult all = priceBermudanLsmcG2pp(b.calibration.model, b.rights);
    const G2ppLsmcResult some = priceBermudanLsmcG2pp(b.calibration.model, {b.rights[1], b.rights[4], b.rights[7]});
    CHECK(all.lower > some.lower + 4.0 * (all.standardError + some.standardError));
    CHECK(priceBermudanLsmcG2pp(b.calibration.model, b.rights).lower == all.lower);
    CHECK_THROWS_AS(priceBermudanLsmcG2pp(b.calibration.model, {b.rights[3], b.rights[1]}), rl::InputError);
    CHECK_THROWS_AS(priceBermudanLsmcG2pp(b.calibration.model, {}), rl::InputError);
}

TEST_CASE("G2++ upper bound: brackets the one-factor PDE with the second factor off, and the gap is small", "[pricers][g2pp][upper]") {
    const CalibratedBermudan one = bermudan10nc1(0.035, true, 0.03);
    const double pde = priceBermudanPde(one.problem).value;
    const CalibratedBermudanG2pp off = calibrateBermudanG2pp(usdMarket(), spotSwap(10, 0.035), 1.0, 0, models::G2pp(0.03, 1.0, 0.01, 0.0, 0.0));
    const G2ppLsmcResult lo = priceBermudanLsmcG2pp(off.calibration.model, off.rights, {.trainingPaths = 20000, .pricingPaths = 40000});
    // The nested estimate of each continuation value is noisy, and max() turns noise into upward bias:
    // with too few inner paths the gap measures the inner simulation, not the rule.
    const G2ppUpperBoundResult coarse = priceBermudanUpperBoundG2pp(off.calibration.model, off.rights, lo, {.outerPaths = 300, .innerPaths = 100});
    const G2ppUpperBoundResult up = priceBermudanUpperBoundG2pp(off.calibration.model, off.rights, lo, {.outerPaths = 300, .innerPaths = 1000});
    CHECK(coarse.gap > up.gap + 2.0 * coarse.gapStandardError);
    INFO("lower " << lo.lower / kBp << " ± " << lo.standardError / kBp << ", upper " << up.upper / kBp << " ± " << up.standardError / kBp << ", PDE " << pde / kBp << " bp");
    CHECK(up.gap >= -3.0 * up.gapStandardError);
    CHECK(up.upper >= pde - 4.0 * up.standardError);
    CHECK(up.gap < 1.0 * kBp);

    // With a genuine second factor the rule regresses on both states and the gap stays small.
    const CalibratedBermudanG2pp two = calibrateBermudanG2pp(usdMarket(), spotSwap(10, 0.035), 1.0, 0, kShape);
    const G2ppLsmcResult lo2 = priceBermudanLsmcG2pp(two.calibration.model, two.rights, {.trainingPaths = 20000, .pricingPaths = 40000});
    const G2ppUpperBoundResult up2 = priceBermudanUpperBoundG2pp(two.calibration.model, two.rights, lo2, {.outerPaths = 300, .innerPaths = 1000});
    CHECK(up2.gap >= -3.0 * up2.gapStandardError);
    CHECK(up2.gap < 1.0 * kBp);
    CHECK_THROWS_AS(priceBermudanUpperBoundG2pp(two.calibration.model, {two.rights[0]}, lo2), rl::InputError);
}
