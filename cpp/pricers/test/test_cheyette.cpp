// The Cheyette model: with no skew it is Hull-White, and must agree with the LGM closed form and PDE;
// with skew it has a smile in the direction of the skew.
#include "rl/pricers/cheyette_mc.hpp"

#include "rl/core/error.hpp"
#include "rl/models/lgm_swaption.hpp"
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
using Catch::Matchers::WithinRel;

namespace {
const double kBp = 1e-4;
const CheyetteOptions kFast{.paths = 40000, .stepsPerYear = 24, .seed = 11};

std::vector<models::ExerciseRight> rights(double rate, bool payer = true) {
    const instruments::Swap swap = spotSwap(10, rate, payer);
    std::vector<models::ExerciseRight> out;
    for (core::Date d : models::coterminalExerciseDates(swap, 1.0, 0, "USD-SIFMA")) out.push_back(models::makeExerciseRight(swap, usdMarket().curves, d));
    return out;
}
} // namespace

TEST_CASE("Cheyette: the deflated bond is a martingale (the y drift term is right)", "[pricers][cheyette]") {
    // An option on a single positive bond is the bond: its price must be D(T), whatever the skew.
    const models::ExerciseRight base = rights(0.035)[4];
    for (double skew : {0.0, 15.0, -10.0}) {
        const models::Cheyette model(0.03, skew, {30.0}, {0.009});
        for (double T : {base.time + 0.5, base.time + 5.0}) {
            models::ExerciseRight bond = base;
            const double df = usdMarket().curves.discount().df(T);
            bond.flows = {{T, 1.0, df}};
            const CheyettePrice p = priceEuropeansCheyette(model, {bond}, kFast)[0];
            CHECK_THAT(p.value, WithinAbs(df, 4.0 * p.standardError + 2e-5 * df));
        }
    }
}

TEST_CASE("Cheyette with no skew is Hull-White: Europeans match the LGM closed form", "[pricers][cheyette]") {
    const double kappa = 0.03, sigma = 0.009;
    const models::Cheyette chey(kappa, 0.0, {30.0}, {sigma});
    const models::Lgm1F lgm = models::Lgm1F::constantSigma(kappa, sigma, 40.0);
    for (double rate : {0.025, 0.035, 0.045}) {
        const auto r = rights(rate);
        const std::vector<CheyettePrice> mc = priceEuropeansCheyette(chey, r, {.paths = 100000, .stepsPerYear = 24, .seed = 5});
        for (std::size_t i = 0; i < r.size(); ++i) {
            const double exact = models::lgmSwaptionPrice(lgm, r[i]);
            INFO("strike " << rate << ", expiry " << r[i].time << ": MC " << mc[i].value / kBp << " bp, closed form " << exact / kBp << " bp");
            CHECK_THAT(mc[i].value, WithinAbs(exact, 4.0 * mc[i].standardError + 0.3 * kBp));
        }
    }
}

TEST_CASE("Cheyette with no skew: the Bermudan matches the LGM PDE", "[pricers][cheyette]") {
    const double kappa = 0.03, sigma = 0.009;
    const auto r = rights(0.035);
    const double pde = priceBermudanPde(BermudanProblem(models::Lgm1F::constantSigma(kappa, sigma, 40.0), r)).value;
    const CheyetteLsmcResult mc = priceBermudanCheyette(models::Cheyette(kappa, 0.0, {30.0}, {sigma}), r, {.trainingPaths = 40000, .pricing = {.paths = 100000, .stepsPerYear = 24, .seed = 3}});
    INFO("LSMC " << mc.lower / kBp << " ± " << mc.standardError / kBp << " bp, PDE " << pde / kBp << " bp");
    CHECK(mc.lower <= pde + 4.0 * mc.standardError + 0.3 * kBp); // a lower bound, up to noise and time-stepping
    CHECK(mc.lower >= pde - 4.0 * mc.standardError - 1.0 * kBp);  // and a good one
}

TEST_CASE("Cheyette calibration reprices the strip, and with no skew finds the Hull-White volatilities", "[pricers][cheyette]") {
    const auto strip = models::makeStrip(rights(0.035), usdMarket().cube);
    const CheyetteCalibration cal = calibrateCheyette(0.03, 0.0, strip, {.paths = 60000, .stepsPerYear = 12, .seed = 9});
    CHECK(cal.ok);
    CHECK(cal.maxAbsErrorBp < 1e-6); // exact on its own paths
    const std::vector<double> hw = models::calibrateLgm(0.03, strip).model.hullWhiteSigmas();
    // Each level absorbs the sampling error of its own swaption's price: a few percent at this path count.
    for (std::size_t i = 0; i < hw.size(); ++i) CHECK_THAT(cal.rows[i].levelBp, WithinRel(hw[i] * 1e4, 0.05));
    // On fresh paths the strip is repriced to its sampling error.
    std::vector<models::ExerciseRight> r;
    for (const auto& ci : strip) r.push_back(ci.right);
    const std::vector<CheyettePrice> fresh = priceEuropeansCheyette(cal.model, r, {.paths = 20000, .stepsPerYear = 12, .seed = 99});
    for (std::size_t i = 0; i < r.size(); ++i) CHECK_THAT(fresh[i].value, WithinAbs(strip[i].marketPrice, 4.0 * fresh[i].standardError + 0.2 * kBp));
    CHECK_THROWS_AS(calibrateCheyette(0.03, 0.0, {}), rl::InputError);
}

TEST_CASE("Cheyette skew produces a smile in its own direction", "[pricers][cheyette]") {
    // Same expiry and swap, three strikes, one simulation; the implied normal vol across strikes.
    const instruments::Swap atmSwap = spotSwap(10, 0.035);
    const core::Date expiry = models::coterminalExerciseDates(atmSwap, 1.0, 0, "USD-SIFMA")[4];
    const double fwd = models::makeExerciseRight(atmSwap, usdMarket().curves, expiry).forwardRate();
    std::vector<models::ExerciseRight> r;
    for (double k : {fwd - 0.01, fwd, fwd + 0.01}) r.push_back(models::makeExerciseRight(spotSwap(10, k), usdMarket().curves, expiry));
    auto smile = [&](double skew) {
        const std::vector<CheyettePrice> p = priceEuropeansCheyette(models::Cheyette(0.03, skew, {30.0}, {0.009}), r, {.paths = 60000, .stepsPerYear = 24, .seed = 4});
        std::vector<double> v;
        for (std::size_t i = 0; i < r.size(); ++i) v.push_back(models::impliedNormalVol(r[i], p[i].value) * 1e4);
        return v;
    };
    const std::vector<double> flat = smile(0.0), up = smile(20.0), down = smile(-15.0);
    INFO("no skew " << flat[0] << " " << flat[1] << " " << flat[2] << "; b>0 " << up[0] << " " << up[1] << " " << up[2] << "; b<0 " << down[0] << " " << down[1] << " " << down[2]);
    CHECK(std::fabs(flat[2] - flat[0]) < 2.0);        // Hull-White: nearly flat in normal vol
    CHECK(up[2] - up[0] > 5.0);                      // vol rising with rates: upward skew
    CHECK(down[2] - down[0] < -5.0);                 // falling with rates: downward skew
}

TEST_CASE("Cheyette Bermudan: at least its best European; one exercise date is the European", "[pricers][cheyette]") {
    const auto r = rights(0.035);
    const models::Cheyette model(0.03, 12.0, {30.0}, {0.009});
    const CheyetteOptions o{.paths = 40000, .stepsPerYear = 24, .seed = 8};
    const std::vector<CheyettePrice> eu = priceEuropeansCheyette(model, r, o);
    double best = 0.0;
    for (const auto& p : eu) best = std::max(best, p.value);
    const CheyetteLsmcResult all = priceBermudanCheyette(model, r, {.trainingPaths = 20000, .pricing = o});
    CHECK(all.lower > best + 10.0 * kBp);
    const CheyetteLsmcResult one = priceBermudanCheyette(model, {r[3]}, {.trainingPaths = 20000, .pricing = o});
    // Same seed and the same time steps to that date: the same paths, so the same number.
    CHECK_THAT(one.lower, WithinAbs(priceEuropeansCheyette(model, {r[3]}, o)[0].value, 1e-12));
    CHECK_THROWS_AS(priceBermudanCheyette(model, {r[2], r[1]}), rl::InputError);
}
