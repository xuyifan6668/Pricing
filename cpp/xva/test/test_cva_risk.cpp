// CVA sensitivities: the adjoint delta against bumping, and the signs a hedge must have.
#include "rl/xva/cva_risk.hpp"

#include "rl/xva/exposure_engine.hpp"
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

TEST_CASE("CS01: negative at the short end, positive at the long end; the ladder adds up to a parallel shift", "[xva][risk]") {
    const SwapExposure e = swapExposureClosedForm(model(), spotSwap(10, 0.035), usdMarket().curves);
    auto cva = [&](const HazardCurve& h) { return creditAdjustment(e.times, e.epe, h, 0.4); };
    const auto ladder = cs01Ladder(cdsQuotes(), 0.4, discount, cva);
    REQUIRE(ladder.size() == 5);
    CHECK(ladder[0].id == "CDS-1Y");
    double sum = 0.0;
    for (const auto& r : ladder) sum += r.value;
    auto all = cdsQuotes();
    for (CdsQuote& q : all) q.spread += 1e-4;
    const double parallel = cva(HazardCurve::bootstrap(all, 0.4, discount)) - cva(counterpartyCurve());
    CHECK(parallel > 0.0);
    CHECK_THAT(sum, WithinRel(parallel, 1e-3));
    // CVA is concave in spreads (default probabilities saturate), so CS01 times the spread level is
    // below the CVA itself, but of its order.
    CHECK(parallel * 125.0 < cva(counterpartyCurve()));
    CHECK(parallel * 125.0 > 0.5 * cva(counterpartyCurve()));
    // Raising one short quote with the others held moves default probability FROM the later years, where
    // the swap's exposure is, INTO the early years, where there is little: CVA falls. A hedge built from
    // one CDS maturity would get the sign of the short end wrong.
    CHECK(ladder[0].value < 0.0);
    CHECK(ladder[1].value < 0.0);
    CHECK(ladder[3].value > 0.0);
    CHECK(ladder[4].value > 0.0);
}

TEST_CASE("Swap CVA delta: one adjoint sweep equals 37 bumps with the model held fixed", "[xva][risk]") {
    const models::UsdMarket& m = usdMarket();
    const HazardCurve hazard = counterpartyCurve();
    for (bool payer : {true, false}) {
        const instruments::Swap swap = spotSwap(10, 0.035, payer);
        const auto adjoint = swapCvaDeltaAdjoint(m, model(), swap, hazard, 0.4);
        // A small central bump: the adjoint is the derivative, a 1 bp one-sided bump is not.
        auto price = [&](const models::UsdMarket& mk) { return swapCvaClosedForm(model(), swap, mk.curves, hazard, 0.4); };
        const auto up = pricers::parDeltaBump(m, price, 0.01), down = pricers::parDeltaBump(m, price, -0.01);
        REQUIRE(adjoint.size() == up.size());
        double total = 0.0;
        for (std::size_t i = 0; i < adjoint.size(); ++i) {
            CHECK(adjoint[i].id == up[i].id);
            CHECK_THAT(adjoint[i].value, WithinAbs(0.5 * (up[i].value + down[i].value), 2e-10));
            total += adjoint[i].value;
        }
        // Higher rates make a payer swap more valuable to its holder: more exposure, more CVA.
        CHECK((payer ? total > 0.0 : total < 0.0));
        CHECK(std::fabs(total) > 1e-7);
    }
}

TEST_CASE("Closed-form swap CVA is the credit adjustment of the closed-form exposure", "[xva][risk]") {
    const instruments::Swap swap = spotSwap(7, 0.03);
    const SwapExposure e = swapExposureClosedForm(model(), swap, usdMarket().curves);
    CHECK_THAT(swapCvaClosedForm(model(), swap, usdMarket().curves, counterpartyCurve(), 0.4), WithinAbs(creditAdjustment(e.times, e.epe, counterpartyCurve(), 0.4), 1e-16));
}
