// The exposure engine, tested against the closed-form swaption strip and on what netting and
// collateral must do whatever the numbers.
#include "rl/xva/exposure_engine.hpp"

#include "rl/core/error.hpp"
#include "rl/models/g2pp_calibration.hpp"
#include "rl/models/g2pp_swaption.hpp"
#include "rl/models/lgm_calibration.hpp"
#include "rl/pricers/exposure.hpp"
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
const market::CurveSet<double>& curves() { return usdMarket().curves; }
/// Mid-period dates of an annual-coupon swap: no coupon rolls off inside a margin period ending here.
std::vector<double> midPeriodTimes(int years) {
    std::vector<double> t{0.0};
    for (int k = 0; k < years; ++k) t.push_back(k + 0.5);
    return t;
}
} // namespace

TEST_CASE("Engine reproduces the closed-form exposure of one swap on its coupon dates", "[xva][engine]") {
    for (bool payer : {true, false}) {
        const instruments::Swap swap = spotSwap(10, 0.035, payer);
        const SwapExposure exact = swapExposureClosedForm(model(), swap, curves());
        const ExposureProfile mc = simulateExposure(model(), {swap}, curves(), exact.times, std::nullopt, {.paths = 40000});
        REQUIRE(mc.times.size() == exact.times.size());
        CHECK_THAT(mc.epe[0], WithinAbs(exact.epe[0], 1e-14));
        for (std::size_t i = 1; i < exact.times.size(); ++i) {
            CHECK_THAT(mc.epe[i], WithinAbs(exact.epe[i], 4.0 * mc.epeStandardError[i]));
            CHECK_THAT(mc.ene[i], WithinAbs(exact.ene[i], 4.0 * mc.epeStandardError[i] + 1e-5));
            CHECK(mc.pfe[i] >= mc.epe[i]);
            CHECK_THAT(mc.grossEpe[i], WithinAbs(mc.epe[i], 1e-14)); // one trade: nothing to net
        }
    }
}

TEST_CASE("Deflated values are martingales: EPE + ENE is today's value of the remaining swap", "[xva][engine]") {
    const instruments::Swap swap = spotSwap(7, 0.030);
    const SwapExposure exact = swapExposureClosedForm(model(), swap, curves());
    const ExposureProfile mc = simulateExposure(model(), {swap}, curves(), exact.times, std::nullopt, {.paths = 40000});
    for (std::size_t i = 1; i < exact.times.size(); ++i) CHECK_THAT(mc.epe[i] + mc.ene[i], WithinAbs(exact.epe[i] + exact.ene[i], 5.0 * mc.epeStandardError[i]));
}

TEST_CASE("Between coupon dates the profile is continuous, and exposure times are what was asked for", "[xva][engine]") {
    const instruments::Swap swap = spotSwap(5, 0.035);
    const std::vector<double> times = exposureTimes({swap}, curves(), 0.25);
    CHECK(times.front() == 0.0);
    CHECK(times.back() < 5.02);
    CHECK(times.size() > 19);
    const ExposureProfile mc = simulateExposure(model(), {swap}, curves(), times, std::nullopt, {.paths = 20000});
    // A coupon in progress is carried at its fixed rate, so no date shows a jump larger than the
    // coupon that has just rolled off plus diffusion.
    for (std::size_t i = 2; i < times.size(); ++i) CHECK(std::fabs(mc.epe[i] - mc.epe[i - 1]) < 0.012);
    CHECK(mc.epe.back() < mc.epe[times.size() / 2]); // the hump: it has come down by the end
}

TEST_CASE("Netting: never worse than gross; offsetting trades cancel; identical trades do not net", "[xva][engine]") {
    const std::vector<double> times = exposureTimes({spotSwap(10, 0.035)}, curves(), 0.5);
    // A payer and the receiver of the same swap: the netted portfolio is worth nothing on every path.
    const ExposureProfile flat = simulateExposure(model(), {spotSwap(10, 0.035, true), spotSwap(10, 0.035, false)}, curves(), times);
    for (std::size_t i = 0; i < times.size(); ++i) {
        CHECK_THAT(flat.epe[i], WithinAbs(0.0, 1e-13));
        if (i > 0) CHECK(flat.grossEpe[i] > 0.005); // without netting both legs would count
    }
    // Two of the same trade: positively dependent, so netting gives nothing.
    const ExposureProfile same = simulateExposure(model(), {spotSwap(10, 0.035), spotSwap(10, 0.035)}, curves(), times);
    for (std::size_t i = 1; i < times.size(); ++i) CHECK_THAT(same.epe[i], WithinAbs(same.grossEpe[i], 1e-12));
    // A mixed book: strictly between.
    const ExposureProfile mixed = simulateExposure(model(), {spotSwap(10, 0.035), spotSwap(5, 0.033, false), spotSwap(7, 0.040)}, curves(), times);
    for (std::size_t i = 1; i < times.size(); ++i) CHECK(mixed.epe[i] <= mixed.grossEpe[i] + 1e-13);
    CHECK(mixed.epe[4] < 0.9 * mixed.grossEpe[4]);
}

TEST_CASE("Collateral: a tighter agreement leaves less exposure; what is left is the margin period", "[xva][engine]") {
    const instruments::Swap swap = spotSwap(10, 0.035);
    const std::vector<double> times = midPeriodTimes(10);
    auto run = [&](const Csa& csa) { return simulateExposure(model(), {swap}, curves(), times, csa); };
    const ExposureProfile perfect = run({.threshold = 0.0, .minimumTransfer = 0.0, .marginPeriodYears = 10.0 / 252.0});
    const ExposureProfile longer = run({.threshold = 0.0, .minimumTransfer = 0.0, .marginPeriodYears = 20.0 / 252.0});
    const ExposureProfile threshold = run({.threshold = 0.01, .minimumTransfer = 0.0, .marginPeriodYears = 10.0 / 252.0});
    const ExposureProfile none = run({.threshold = 1e9, .minimumTransfer = 0.0, .marginPeriodYears = 10.0 / 252.0});
    for (std::size_t i = 1; i < times.size(); ++i) {
        CHECK(perfect.collateralisedEpe[i] > 0.0);                            // never zero: the margin period
        CHECK(perfect.collateralisedEpe[i] < 0.5 * perfect.epe[i]);
        CHECK(longer.collateralisedEpe[i] > perfect.collateralisedEpe[i]);   // grows with the margin period
        CHECK(threshold.collateralisedEpe[i] > perfect.collateralisedEpe[i]); // and with the threshold
        CHECK(threshold.collateralisedEpe[i] < threshold.epe[i]);
        CHECK_THAT(none.collateralisedEpe[i], WithinAbs(none.epe[i], 1e-14)); // an unreachable threshold is no agreement
        // Over a short period the residual exposure scales like the square root of its length.
        CHECK_THAT(longer.collateralisedEpe[i] / perfect.collateralisedEpe[i], WithinAbs(std::sqrt(2.0), 0.08));
    }
    // The residual exposure follows the swap's remaining sensitivity down to nothing at maturity.
    CHECK(perfect.collateralisedEpe.back() < 0.2 * perfect.collateralisedEpe[1]);
    // The uncollateralised numbers do not depend on the agreement.
    CHECK_THAT(perfect.epe[5], WithinAbs(threshold.epe[5], 1e-14));
}

TEST_CASE("Collateral: a coupon rolling off inside the margin period leaves a spike", "[xva][engine]") {
    // On a coupon date the portfolio's value jumps by the net coupon, but the collateral held was
    // called on the value before the jump. Whoever was owed the coupon is left exposed for it.
    const instruments::Swap swap = spotSwap(10, 0.035);
    const std::vector<double> starts = periodStartTimes(swap, curves()); // 0, spot, then each anniversary
    std::vector<double> times{0.0};
    for (std::size_t k = 2; k < starts.size(); ++k) {
        times.push_back(starts[k] - 0.5); // mid-period
        times.push_back(starts[k]);       // just after a coupon has rolled off
    }
    const ExposureProfile p = simulateExposure(model(), {swap}, curves(), times, Csa{0.0, 0.0, 10.0 / 252.0});
    // The first coupon's rate is known today, so its roll-off is no surprise; every later one is.
    for (std::size_t i = 3; i + 1 < times.size(); i += 2) {
        INFO("t = " << times[i + 1] << ": " << p.collateralisedEpe[i] << " mid-period, " << p.collateralisedEpe[i + 1] << " after the coupon");
        CHECK(p.collateralisedEpe[i + 1] > 1.3 * p.collateralisedEpe[i]);
    }
}

TEST_CASE("Minimum transfer amount changes the residual exposure only at second order", "[xva][engine]") {
    // A call that is not made leaves the portfolio's small value uncollateralised in both directions,
    // so to first order it adds as much as it removes.
    const instruments::Swap swap = spotSwap(10, 0.035);
    const std::vector<double> times = midPeriodTimes(10);
    const ExposureProfile without = simulateExposure(model(), {swap}, curves(), times, Csa{0.0, 0.0, 10.0 / 252.0});
    const ExposureProfile with = simulateExposure(model(), {swap}, curves(), times, Csa{0.0, 0.005, 10.0 / 252.0});
    for (std::size_t i = 1; i < times.size(); ++i) {
        CHECK(with.collateralisedEpe[i] != without.collateralisedEpe[i]);
        CHECK_THAT(with.collateralisedEpe[i], WithinAbs(without.collateralisedEpe[i], 0.25 * without.collateralisedEpe[i] + 2e-4));
    }
}

TEST_CASE("Engine rejects inputs it cannot value", "[xva][engine]") {
    const instruments::Swap swap = spotSwap(5, 0.035);
    CHECK_THROWS_AS(simulateExposure(model(), {}, curves(), {0.0, 1.0}), rl::InputError);
    CHECK_THROWS_AS(simulateExposure(model(), {swap}, curves(), {1.0, 2.0}), rl::InputError);
    CHECK_THROWS_AS(simulateExposure(model(), {swap}, curves(), {0.0, 2.0, 1.0}), rl::InputError);
    CHECK_THROWS_AS(simulateExposure(model(), {swap}, curves(), {0.0, 1.0}, Csa{-1.0, 0.0, 0.04}), rl::InputError);
    const instruments::Swap ff = instruments::OisSwapBuilder::usdFedFunds(usdMarket().asOf).id("ff").tenor(core::Tenor::years(5)).fixedRate(0.035).build();
    CHECK_THROWS_AS(simulateExposure(model(), {ff}, curves(), {0.0, 1.0}), rl::InputError);
}

namespace {
/// A two-factor shape with its volatility scale calibrated to the 10-year co-terminal strip at 3.50%.
models::G2pp twoFactor(double rho, double sigma2 = 0.009, double kappa1 = 0.40) {
    const instruments::Swap swap = spotSwap(10, 0.035);
    std::vector<models::ExerciseRight> rights;
    for (core::Date d : models::coterminalExerciseDates(swap, 1.0, 0, "USD-SIFMA")) rights.push_back(models::makeExerciseRight(swap, curves(), d));
    return models::calibrateG2ppToStrip(models::G2pp(kappa1, 0.03, 0.012, sigma2, rho), models::makeStrip(rights, usdMarket().cube)).model;
}
} // namespace

TEST_CASE("Two-factor engine reproduces the two-factor swaption strip on coupon dates", "[xva][engine][g2pp]") {
    const models::G2pp g = twoFactor(-0.7);
    for (bool payer : {true, false}) {
        const instruments::Swap swap = spotSwap(7, 0.034, payer);
        const std::vector<double> times = periodStartTimes(swap, curves());
        const ExposureProfile mc = simulateExposure(g, {swap}, curves(), times, std::nullopt, {.paths = 40000});
        std::size_t i = 1;
        double last = 0.0;
        for (const auto& cf : swap.fixedLeg().cashflows()) {
            const core::Date start = std::get<instruments::FixedCoupon>(cf).accrualStart;
            const double t = curves().discount().time(start);
            if (t <= last + 1e-9) continue;
            last = t;
            const std::vector<models::BondFlow> flows = models::swapAsBonds(swap, curves(), start);
            double forward = 0.0;
            for (const auto& f : flows) forward += f.amount * f.df;
            const double exact = models::g2ppOptionValue(g, t, flows);
            CHECK_THAT(mc.epe[i], WithinAbs(exact, 4.0 * mc.epeStandardError[i]));
            CHECK_THAT(mc.epe[i] + mc.ene[i], WithinAbs(forward, 5.0 * mc.epeStandardError[i])); // deflated values are martingales
            ++i;
        }
        CHECK(i == times.size());
    }
}

TEST_CASE("Two-factor engine with the second factor off is the one-factor engine", "[xva][engine][g2pp]") {
    const instruments::Swap swap = spotSwap(10, 0.035);
    const SwapExposure exact = swapExposureClosedForm(model(), swap, curves());
    // Same mean reversion as the one-factor reference and no second factor: the same model.
    const ExposureProfile mc = simulateExposure(twoFactor(0.0, 0.0, 0.03), {swap}, curves(), exact.times, std::nullopt, {.paths = 40000});
    for (std::size_t i = 1; i < exact.times.size(); ++i) CHECK_THAT(mc.epe[i], WithinAbs(exact.epe[i], 4.0 * mc.epeStandardError[i]));
}

TEST_CASE("Factor correlation changes how trades of different maturities net, not the calibrated swap", "[xva][engine][g2pp]") {
    const instruments::Swap ten = spotSwap(10, 0.035);
    const std::vector<double> couponDates = periodStartTimes(ten, curves()), midDates = midPeriodTimes(10);
    const std::vector<instruments::Swap> steepener{spotSwap(10, 0.035, true, 1.0), spotSwap(2, 0.035, false, 4.5)}; // long-end payer against a short-end receiver
    std::vector<double> tenLow, tenHigh;
    double bookLow = 0.0, bookHigh = 0.0;
    for (double rho : {-0.9, 0.5}) {
        const models::G2pp g = twoFactor(rho);
        const ExposureProfile a = simulateExposure(g, {ten}, curves(), couponDates, std::nullopt, {.paths = 40000}), b = simulateExposure(g, steepener, curves(), midDates);
        (rho < 0.0 ? tenLow : tenHigh) = a.epe;
        (rho < 0.0 ? bookLow : bookHigh) = b.epe[1];
        for (std::size_t i = 1; i < midDates.size(); ++i) CHECK(b.epe[i] <= b.grossEpe[i] + 1e-13);
    }
    // The 10-year swap's exposure on its coupon dates IS the calibration set, so it does not move ...
    for (std::size_t i = 2; i < couponDates.size(); ++i) CHECK_THAT(tenLow[i], WithinAbs(tenHigh[i], 0.02 * tenHigh[i]));
    // ... while a position in the slope of the curve is exposed to exactly what the correlation controls.
    CHECK(bookLow > 1.25 * bookHigh);
}

TEST_CASE("Two-factor engine with collateral: residual exposure is positive and small", "[xva][engine][g2pp]") {
    const std::vector<double> times = midPeriodTimes(10);
    const ExposureProfile p = simulateExposure(twoFactor(-0.7), {spotSwap(10, 0.035)}, curves(), times, Csa{0.0, 0.0, 10.0 / 252.0});
    for (std::size_t i = 1; i < times.size(); ++i) {
        CHECK(p.collateralisedEpe[i] > 0.0);
        CHECK(p.collateralisedEpe[i] < 0.5 * p.epe[i]);
    }
}

TEST_CASE("PFE: the engine's quantile matches the closed form for one swap", "[xva][engine][pfe]") {
    // A swap's value is monotone in the state, so the 97.5% quantile of its positive exposure is its
    // value at the matching quantile of the state: a closed form the simulated quantile must reproduce.
    for (bool payer : {true, false}) {
        const instruments::Swap swap = spotSwap(10, 0.035, payer);
        const std::vector<pricers::ExposurePoint> exact = pricers::swapExposureProfile(model(), swap, curves());
        const std::vector<double> times = periodStartTimes(swap, curves());
        const ExposureProfile mc = simulateExposure(model(), {swap}, curves(), times, std::nullopt, {.paths = 40000, .quantile = 0.975});
        REQUIRE(times.size() == exact.size() + 1);
        for (std::size_t i = 2; i < times.size(); ++i) {
            INFO("t = " << times[i] << ": engine " << mc.pfe[i] / 1e-4 << " bp, closed form " << exact[i - 1].potentialExposure / 1e-4 << " bp");
            CHECK_THAT(times[i], WithinAbs(exact[i - 1].time, 1e-9));
            // The sampling error of a 97.5% quantile from 40,000 draws is about 1.5% of its level here.
            CHECK_THAT(mc.pfe[i], WithinRel(exact[i - 1].potentialExposure, 0.04));
        }
    }
}

TEST_CASE("Monte Carlo CVA: unbiased across seeds, and its error falls like one over root paths", "[xva][engine][convergence]") {
    const instruments::Swap swap = spotSwap(10, 0.035);
    const SwapExposure exact = swapExposureClosedForm(model(), swap, curves());
    const HazardCurve hazard = counterpartyCurve();
    const double truth = creditAdjustment(exact.times, exact.epe, hazard, 0.4);
    auto spread = [&](std::size_t paths, double& mean) {
        std::vector<double> cva;
        for (std::uint64_t seed = 1; seed <= 12; ++seed) {
            const ExposureProfile p = simulateExposure(model(), {swap}, curves(), exact.times, std::nullopt, {.paths = paths, .seed = seed});
            cva.push_back(creditAdjustment(p.times, p.epe, hazard, 0.4));
        }
        mean = 0.0;
        for (double v : cva) mean += v;
        mean /= static_cast<double>(cva.size());
        double var = 0.0;
        for (double v : cva) var += (v - mean) * (v - mean);
        return std::sqrt(var / static_cast<double>(cva.size() - 1));
    };
    double meanSmall = 0.0, meanLarge = 0.0;
    const double sdSmall = spread(1000, meanSmall), sdLarge = spread(4000, meanLarge);
    INFO("closed form " << truth / 1e-4 << " bp; 1,000 paths " << meanSmall / 1e-4 << " ± " << sdSmall / 1e-4 << "; 4,000 paths " << meanLarge / 1e-4 << " ± " << sdLarge / 1e-4);
    // Unbiased: the seed average is within its own error of the closed form.
    CHECK_THAT(meanSmall, WithinAbs(truth, 4.0 * sdSmall / std::sqrt(12.0)));
    CHECK_THAT(meanLarge, WithinAbs(truth, 4.0 * sdLarge / std::sqrt(12.0)));
    // Four times the paths: half the error (with twelve seeds the ratio itself is only known to about 30%).
    CHECK_THAT(sdSmall / sdLarge, WithinAbs(2.0, 0.8));
}
