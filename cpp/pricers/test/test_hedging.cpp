// The delta hedge and the P&L explain, tested on what theory says they must do.
#include "rl/pricers/hedging.hpp"

#include "rl/core/error.hpp"
#include "rl/pricers/bermudan_pde.hpp"
#include "test_support.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>

using namespace rl;
using namespace rl::pricers;
using namespace rl::pricers::test;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace {
const double kBp = 1e-4;
}

TEST_CASE("Delta hedge: the model's own price and the PDE agree; the hedge removes most of the risk", "[pricers][hedging]") {
    const CalibratedBermudan b = bermudan10nc1();
    const HedgingResult naked = simulateDeltaHedge(b.problem.model(), b.problem.rights(), {.paths = 4000, .hedge = false});
    const HedgingResult hedged = simulateDeltaHedge(b.problem.model(), b.problem.rights(), {.paths = 4000, .rebalancesPerYear = 52});
    PdeOptions fine;
    fine.spaceSteps = fine.timeSteps = 1600; // the default grid is 0.08 bp from converged here
    CHECK_THAT(hedged.truePrice, WithinAbs(priceBermudanPde(b.problem, fine).value, 0.02 * kBp));
    CHECK(hedged.premium == hedged.truePrice); // vol scale 1: the hedger's model is the truth
    CHECK(hedged.pnlStdDev < 0.1 * naked.pnlStdDev);
    // With the right model the hedged P&L has mean zero.
    CHECK_THAT(hedged.meanPnl, WithinAbs(0.0, 4.0 * hedged.meanPnlStandardError));
}

TEST_CASE("Delta hedge: the error falls like the square root of the rebalancing interval", "[pricers][hedging]") {
    const CalibratedBermudan b = bermudan10nc1();
    const HedgingResult monthly = simulateDeltaHedge(b.problem.model(), b.problem.rights(), {.paths = 3000, .rebalancesPerYear = 12});
    const HedgingResult weekly = simulateDeltaHedge(b.problem.model(), b.problem.rights(), {.paths = 3000, .rebalancesPerYear = 48});
    // Four times as often: half the error.
    CHECK_THAT(monthly.pnlStdDev / weekly.pnlStdDev, WithinAbs(2.0, 0.3));
}

TEST_CASE("Delta hedge with the wrong vol: the mean P&L is the mispricing", "[pricers][hedging]") {
    const CalibratedBermudan b = bermudan10nc1();
    for (double s : {0.8, 1.2}) {
        const HedgingResult r = simulateDeltaHedge(b.problem.model(), b.problem.rights(), {.paths = 4000, .rebalancesPerYear = 52, .volScale = s});
        CHECK((s < 1.0 ? r.premium < r.truePrice : r.premium > r.truePrice));
        // Sold too cheap loses on average, too dear gains; the hedge does not undo the price.
        CHECK_THAT(r.meanPnl, WithinAbs(r.premium - r.truePrice, 4.0 * r.meanPnlStandardError + 0.05 * std::fabs(r.premium - r.truePrice)));
    }
    CHECK_THROWS_AS(simulateDeltaHedge(b.problem.model(), {}), rl::InputError);
}

TEST_CASE("P&L explain: delta leaves a second-order residual, delta plus gamma a third-order one", "[pricers][hedging][pnl]") {
    const models::UsdMarket& m = usdMarket();
    const instruments::Swap trade = spotSwap(10, 0.035);
    PdeOptions fine;
    fine.spaceSteps = fine.timeSteps = 1600; // a fine grid: second-order risk is a difference of differences
    auto price = [&](const models::UsdMarket& mk) { return priceBermudanPde(calibrateBermudan(mk, trade, 1.0, 0, 0.03).problem, fine).value; };
    const Sensitivities s = computeSensitivities(m, price);
    CHECK(s.parallelGamma > 0.0); // a long option is long gamma
    CHECK(s.vega > 0.0);
    auto residuals = [&](double h) {
        MarketMove move;
        move.quoteShiftBp.assign(m.problem.size(), h);
        const PnlExplain e = explainPnl(m, s, price, move);
        return std::pair{e.actual - e.delta, e.unexplained()};
    };
    const auto [d10, g10] = residuals(10.0);
    const auto [d20, g20] = residuals(20.0);
    INFO("delta-only residual " << d10 / kBp << ", " << d20 / kBp << " bp; with gamma " << g10 / kBp << ", " << g20 / kBp << " bp");
    CHECK_THAT(d20 / d10, WithinAbs(4.0, 0.8));   // second order
    CHECK(std::fabs(g10) < 0.15 * std::fabs(d10)); // gamma explains most of what delta leaves
    CHECK(std::fabs(g20) < 0.15 * std::fabs(d20));
    // A pure vol move is explained by vega to second order.
    MarketMove vol;
    vol.volShiftBp = 2.0;
    const PnlExplain ev = explainPnl(m, s, price, vol);
    CHECK(std::fabs(ev.unexplained()) < 0.05 * std::fabs(ev.actual));
}
