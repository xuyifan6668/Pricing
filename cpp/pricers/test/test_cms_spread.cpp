// The CMS spread option: limits where the answer is known, and what correlation must do to it.
#include "rl/pricers/cms_spread.hpp"

#include "rl/core/error.hpp"
#include "rl/math/rng.hpp"
#include "rl/models/g2pp_calibration.hpp"
#include "rl/models/g2pp_swaption.hpp"
#include "rl/models/lgm_swaption.hpp"
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
/// The 5-year-forward 2-year and 10-year ATM swaptions, and D(5y).
struct Setup {
    std::vector<models::CalibrationInstrument> ins;
    double df;
};
const Setup& setup() {
    static const Setup s = [] {
        const models::UsdMarket& m = usdMarket();
        Setup out{models::makeAtmMatrix(m.asOf, m.curves, m.cube, {5}, {2, 10}), 0.0};
        out.df = m.curves.discount().df(out.ins[0].right.time);
        return out;
    }();
    return s;
}
const models::ExerciseRight& shortSwap() { return setup().ins[0].right; }
const models::ExerciseRight& longSwap() { return setup().ins[1].right; }
const models::G2pp kModel(0.40, 0.03, 0.012, 0.009, -0.7);
} // namespace

TEST_CASE("CMS spread: the model's rate volatilities are the swaption vols it prices", "[pricers][cms]") {
    const CmsSpreadResult r = priceCmsSpread(kModel, longSwap(), shortSwap(), 0.0, true, setup().df);
    // The standard deviation of each rate at expiry is, to convexity, the normal vol of its swaption.
    CHECK_THAT(r.longVolBp, WithinRel(models::impliedNormalVol(longSwap(), models::g2ppSwaptionPrice(kModel, longSwap())) * 1e4, 0.01));
    CHECK_THAT(r.shortVolBp, WithinRel(models::impliedNormalVol(shortSwap(), models::g2ppSwaptionPrice(kModel, shortSwap())) * 1e4, 0.01));
    // The expected rates are the forwards plus a small convexity adjustment.
    CHECK(r.longRate > longSwap().forwardRate());
    CHECK(r.shortRate > shortSwap().forwardRate());
    CHECK_THAT(r.longRate, WithinAbs(longSwap().forwardRate(), 20.0 * kBp));
    CHECK_THAT(r.shortRate, WithinAbs(shortSwap().forwardRate(), 20.0 * kBp));
    CHECK(r.rateCorrelation > 0.3);
    CHECK(r.rateCorrelation < 0.98);
    CHECK_THAT(r.spreadVolBp * r.spreadVolBp, WithinRel(r.longVolBp * r.longVolBp + r.shortVolBp * r.shortVolBp - 2.0 * r.rateCorrelation * r.longVolBp * r.shortVolBp, 1e-9));
}

TEST_CASE("CMS spread: parity, intrinsic limits and convergence of the integral", "[pricers][cms]") {
    const double df = setup().df;
    const CmsSpreadResult atm = priceCmsSpread(kModel, longSwap(), shortSwap(), 0.0, true, df);
    for (double k : {-0.005, 0.0, atm.forwardSpread, 0.01}) {
        const CmsSpreadResult c = priceCmsSpread(kModel, longSwap(), shortSwap(), k, true, df), p = priceCmsSpread(kModel, longSwap(), shortSwap(), k, false, df);
        CHECK_THAT(c.price - p.price, WithinAbs(df * (c.forwardSpread - k), 1e-10)); // call − put = the forward
        CHECK(c.price >= df * std::max(c.forwardSpread - k, 0.0) - 1e-10);
    }
    // Far in and far out of the money.
    CHECK_THAT(priceCmsSpread(kModel, longSwap(), shortSwap(), -0.20, true, df).price, WithinAbs(df * (atm.forwardSpread + 0.20), 1e-9));
    CHECK_THAT(priceCmsSpread(kModel, longSwap(), shortSwap(), 0.20, true, df).price, WithinAbs(0.0, 1e-12));
    // A spread of a rate with itself is zero in every state.
    const CmsSpreadResult self = priceCmsSpread(kModel, longSwap(), longSwap(), -0.01, true, df);
    CHECK_THAT(self.price, WithinAbs(df * 0.01, 1e-10));
    CHECK_THAT(self.spreadVolBp, WithinAbs(0.0, 1e-4));
    CHECK_THAT(priceCmsSpread(kModel, longSwap(), shortSwap(), atm.forwardSpread, true, df, 801).price, WithinAbs(priceCmsSpread(kModel, longSwap(), shortSwap(), atm.forwardSpread, true, df).price, 0.002 * kBp));
    CHECK_THROWS_AS(priceCmsSpread(kModel, longSwap(), shortSwap(), 0.0, true, df, 400), rl::InputError);
    const models::ExerciseRight other = models::makeAtmMatrix(usdMarket().asOf, usdMarket().curves, usdMarket().cube, {3}, {2})[0].right;
    CHECK_THROWS_AS(priceCmsSpread(kModel, longSwap(), other, 0.0, true, df), rl::InputError);
}

TEST_CASE("CMS spread: agrees with a Monte Carlo of the same payoff", "[pricers][cms]") {
    const double df = setup().df, t = longSwap().time;
    const CmsSpreadResult r = priceCmsSpread(kModel, longSwap(), shortSwap(), 0.002, true, df);
    // Draw the state at expiry and evaluate both par rates from deflated bonds directly.
    const models::G2pp::Covariance c = kModel.covariance(t);
    const double sd1 = std::sqrt(c.c11), beta = c.c12 / c.c11, sd2 = std::sqrt(c.c22 - c.c12 * beta);
    auto rate = [&](const models::ExerciseRight& right, double x1, double x2) {
        auto sum = [&](const std::vector<models::BondFlow>& flows) {
            double v = 0.0;
            for (const auto& f : flows) {
                const auto l = kModel.loading(f.time, t);
                v += f.amount * f.df * std::exp(-l.h1 * x1 - l.h2 * x2 - 0.5 * models::G2pp::quadraticForm(l, c));
            }
            return v;
        };
        return right.strike + sum(right.flows) / sum(right.annuity);
    };
    math::Rng rng(11);
    const int n = 200000;
    double sum = 0.0, sumSq = 0.0;
    for (int i = 0; i < n; ++i) {
        const double x1 = sd1 * rng.gaussian(), x2 = beta * x1 + sd2 * rng.gaussian();
        const double v = df * std::max(rate(longSwap(), x1, x2) - rate(shortSwap(), x1, x2) - 0.002, 0.0);
        sum += v;
        sumSq += v * v;
    }
    const double mean = sum / n, se = std::sqrt((sumSq / n - mean * mean) / n);
    CHECK_THAT(r.price, WithinAbs(mean, 4.0 * se));
}

TEST_CASE("CMS spread: with both swaption vols matched, the price is a reading of the correlation", "[pricers][cms]") {
    const double df = setup().df;
    const auto& ins = setup().ins;
    double previousPrice = 1e9, previousCorrelation = -2.0;
    for (double rho : {-0.95, -0.85, -0.75}) {
        // Two volatility ratios reprice the pair at each rho; follow the branch with the smaller second factor.
        const models::G2pp model = models::calibrateG2ppToPair(0.40, 0.03, rho, ins[0], ins[1], 0.4);
        const CmsSpreadResult atm = priceCmsSpread(model, longSwap(), shortSwap(), longSwap().forwardRate() - shortSwap().forwardRate(), true, df);
        CHECK(atm.rateCorrelation > previousCorrelation); // more factor correlation, more rate correlation
        CHECK(atm.price < previousPrice);                 // and a cheaper option on their difference
        previousPrice = atm.price;
        previousCorrelation = atm.rateCorrelation;
    }
    // One factor, matched to the same two vols: the rates move together and little is left.
    const models::G2pp one = models::calibrateOneFactorToPair(ins[0], ins[1]);
    const CmsSpreadResult r = priceCmsSpread(one, longSwap(), shortSwap(), longSwap().forwardRate() - shortSwap().forwardRate(), true, df);
    CHECK(r.rateCorrelation > 0.999);
    CHECK(r.spreadVolBp < 2.0);
    CHECK(r.price < 0.5 * previousPrice);
}
