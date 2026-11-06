// CVA, DVA and FVA: closed form against simulation, and the signs and limits each must have.
#include "rl/xva/xva.hpp"

#include "rl/core/error.hpp"
#include "rl/xva/exposure_engine.hpp"
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
const market::CurveSet<double>& curves() { return usdMarket().curves; }
} // namespace

TEST_CASE("Credit adjustment: constant exposure gives loss-given-default times default probability", "[xva]") {
    const HazardCurve h = HazardCurve::flat(0.02);
    const std::vector<double> times{0.0, 1.0, 2.5, 5.0}, flat(4, 0.03);
    CHECK_THAT(creditAdjustment(times, flat, h, 0.4), WithinRel(0.6 * 0.03 * (1.0 - std::exp(-0.1)), 1e-13));
    CHECK(creditAdjustment(times, flat, HazardCurve::flat(0.0), 0.4) == 0.0);
    CHECK(creditAdjustment(times, flat, h, 1.0) == 0.0); // full recovery: nothing to lose
    CHECK_THAT(fundingAdjustment(times, flat, 0.005), WithinRel(0.005 * 0.03 * 5.0, 1e-13));
    CHECK_THROWS_AS(creditAdjustment({0.0, 1.0}, {0.1}, h, 0.4), rl::InputError);
    CHECK_THROWS_AS(creditAdjustment({0.0, 1.0}, {0.1, 0.1}, h, 1.5), rl::InputError);
}

TEST_CASE("Swap CVA: Monte Carlo agrees with the closed-form swaption strip", "[xva]") {
    const HazardCurve hazard = counterpartyCurve();
    for (bool payer : {true, false}) {
        const instruments::Swap swap = spotSwap(10, 0.035, payer);
        const SwapExposure exact = swapExposureClosedForm(model(), swap, curves());
        const ExposureProfile mc = simulateExposure(model(), {swap}, curves(), exact.times, std::nullopt, {.paths = 40000});
        const double closed = creditAdjustment(exact.times, exact.epe, hazard, 0.4), simulated = creditAdjustment(mc.times, mc.epe, hazard, 0.4);
        CHECK(closed > 5e-4);  // a few basis points of notional on a 10-year swap at these spreads
        CHECK(closed < 5e-3);
        CHECK_THAT(simulated, WithinRel(closed, 0.02));
    }
}

TEST_CASE("One side's CVA exposure is the other side's DVA exposure", "[xva]") {
    const SwapExposure pay = swapExposureClosedForm(model(), spotSwap(10, 0.035, true), curves());
    const SwapExposure rec = swapExposureClosedForm(model(), spotSwap(10, 0.035, false), curves());
    for (std::size_t i = 0; i < pay.times.size(); ++i) CHECK_THAT(pay.epe[i], WithinAbs(-rec.ene[i], 1e-13));
    XvaInputs in;
    in.counterparty = counterpartyCurve();
    in.own = counterpartyCurve();
    in.fundingSpread = 0.006;
    const XvaResult a = computeXva(pay.times, pay.epe, pay.ene, in), b = computeXva(rec.times, rec.epe, rec.ene, in);
    CHECK_THAT(a.cva, WithinRel(b.dva, 1e-12));
    CHECK_THAT(a.dva, WithinRel(b.cva, 1e-12));
    CHECK_THAT(a.fca, WithinRel(b.fba, 1e-12));
    CHECK_THAT(a.fva(), WithinRel(-b.fva(), 1e-12));
    CHECK(a.cva > 0.0);
    CHECK(a.dva > 0.0);
}

TEST_CASE("CVA rises with spreads, falls with recovery at fixed hazard, and vanishes with collateral", "[xva]") {
    const instruments::Swap swap = spotSwap(10, 0.035);
    const SwapExposure e = swapExposureClosedForm(model(), swap, curves());
    const double base = creditAdjustment(e.times, e.epe, counterpartyCurve(), 0.4);
    auto wide = cdsQuotes();
    for (CdsQuote& q : wide) q.spread *= 2.0;
    const double doubled = creditAdjustment(e.times, e.epe, HazardCurve::bootstrap(wide, 0.4, discount), 0.4);
    CHECK(doubled > 1.8 * base);
    CHECK(doubled < 2.0 * base); // default probabilities are concave in the spread
    CHECK(creditAdjustment(e.times, e.epe, counterpartyCurve(), 0.6) < base);

    std::vector<double> times{0.0};
    for (int k = 0; k < 10; ++k) times.push_back(k + 0.5);
    const ExposureProfile p = simulateExposure(model(), {swap}, curves(), times, Csa{0.0, 0.0, 10.0 / 252.0});
    const double open = creditAdjustment(p.times, p.epe, counterpartyCurve(), 0.4);
    const double collateralised = creditAdjustment(p.times, p.collateralisedEpe, counterpartyCurve(), 0.4);
    CHECK(collateralised < 0.35 * open);
    CHECK(collateralised > 0.0);
}
