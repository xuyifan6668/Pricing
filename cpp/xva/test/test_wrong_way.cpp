// Wrong-way risk: the independent limit, the sign of the effect, and the curve it must still reprice.
#include "rl/xva/wrong_way.hpp"

#include "rl/core/error.hpp"
#include "rl/xva/xva.hpp"
#include "test_support.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>

using namespace rl;
using namespace rl::xva;
using namespace rl::xva::test;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace {
ExposureProfile paths(bool payer) {
    const instruments::Swap swap = spotSwap(10, 0.035, payer);
    return simulateExposure(model(), {swap}, usdMarket().curves, exposureTimes({swap}, usdMarket().curves, 0.5), std::nullopt, {.paths = 40000, .keepPaths = true});
}
} // namespace

TEST_CASE("Wrong-way risk: no dependence is the independent CVA", "[xva][wwr]") {
    const ExposureProfile p = paths(true);
    const WrongWayResult r = wrongWayCva(p, model(), counterpartyCurve(), 0.4, 0.0);
    // The same paths, the same trapezoid, the curve's own default probabilities.
    CHECK_THAT(r.cva, WithinRel(r.independentCva, 1e-10));
    CHECK_THAT(r.independentCva, WithinRel(creditAdjustment(p.times, p.epe, counterpartyCurve(), 0.4), 1e-14));
    // With no dependence, exposure given default is just exposure, and the hazard levels are the curve's
    // forward hazard rates.
    for (std::size_t j = 2; j < r.times.size(); ++j) CHECK_THAT(r.epeGivenDefault[j], WithinRel(r.epe[j], 1e-10));
    CHECK_THAT(r.hazardLevel[1], WithinRel(counterpartyCurve().hazards()[0], 1e-9));
    CHECK_THAT(r.hazardLevel.back(), WithinRel(counterpartyCurve().hazards().back(), 1e-9));
}

TEST_CASE("Wrong-way risk: the CDS curve is repriced whatever the dependence", "[xva][wwr]") {
    const ExposureProfile p = paths(true);
    for (double beta : {-1.0, -0.5, 0.5, 1.0}) {
        const WrongWayResult r = wrongWayCva(p, model(), counterpartyCurve(), 0.4, beta);
        for (std::size_t j = 1; j < r.times.size(); ++j) {
            CHECK_THAT(r.survivalError[j], WithinAbs(0.0, 1e-12));
            CHECK(r.hazardLevel[j] > 0.0);
        }
    }
}

TEST_CASE("Wrong-way risk: credit that worsens with rates hurts the payer and helps the receiver", "[xva][wwr]") {
    const ExposureProfile payer = paths(true), receiver = paths(false);
    const HazardCurve hazard = counterpartyCurve();
    double previousPayer = 0.0, previousReceiver = 1e9;
    for (double beta : {-1.0, -0.5, 0.0, 0.5, 1.0}) {
        const WrongWayResult a = wrongWayCva(payer, model(), hazard, 0.4, beta), b = wrongWayCva(receiver, model(), hazard, 0.4, beta);
        CHECK(a.cva > previousPayer);       // the payer is owed more when rates are high
        CHECK(b.cva < previousReceiver);    // the receiver when they are low
        previousPayer = a.cva;
        previousReceiver = b.cva;
        if (beta > 0.0) {
            CHECK(a.cva > a.independentCva + 4.0 * a.standardError);
            CHECK(b.cva < b.independentCva);
            // Exposure conditional on default is above the unconditional one for the payer, below for the receiver.
            CHECK(a.epeGivenDefault[10] > a.epe[10]);
            CHECK(b.epeGivenDefault[10] < b.epe[10]);
        }
    }
    // A payer facing beta is a receiver facing -beta: the two effects are of the same order.
    const double up = wrongWayCva(payer, model(), hazard, 0.4, 1.0).cva / wrongWayCva(payer, model(), hazard, 0.4, 0.0).cva;
    const double down = wrongWayCva(receiver, model(), hazard, 0.4, -1.0).cva / wrongWayCva(receiver, model(), hazard, 0.4, 0.0).cva;
    CHECK(up > 1.2);
    CHECK(down > 1.2);
}

TEST_CASE("Wrong-way risk rejects profiles it cannot use", "[xva][wwr]") {
    const instruments::Swap swap = spotSwap(5, 0.035);
    const ExposureProfile noPaths = simulateExposure(model(), {swap}, usdMarket().curves, {0.0, 1.0, 2.0});
    CHECK_THROWS_AS(wrongWayCva(noPaths, model(), counterpartyCurve(), 0.4, 0.5), rl::InputError);
    CHECK_THROWS_AS(wrongWayCva(paths(true), model(), counterpartyCurve(), 1.5, 0.5), rl::InputError);
}
