// The two-factor model and European swaptions under it, tested against limits and independent
// formulas rather than stored numbers.
#include "rl/models/g2pp.hpp"
#include "rl/models/g2pp_swaption.hpp"

#include "rl/core/error.hpp"
#include "rl/math/normal.hpp"
#include "rl/math/quadrature.hpp"
#include "rl/models/lgm.hpp"
#include "rl/models/lgm_swaption.hpp"
#include "test_support.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>

using namespace rl::models;
using namespace rl::models::test;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace {
const G2pp kModel(0.35, 0.02, 0.011, 0.009, -0.75); // a fast and a slow factor, negatively correlated
}

TEST_CASE("G2++ covariance: closed form for constant parameters, additive across scale intervals", "[models][g2pp]") {
    const double t = 3.0;
    const G2pp::Covariance c = kModel.covariance(t);
    CHECK_THAT(c.c11, WithinRel(0.011 * 0.011 * std::expm1(2 * 0.35 * t) / (2 * 0.35), 1e-13));
    CHECK_THAT(c.c22, WithinRel(0.009 * 0.009 * std::expm1(2 * 0.02 * t) / (2 * 0.02), 1e-13));
    CHECK_THAT(c.c12, WithinRel(-0.75 * 0.011 * 0.009 * std::expm1(0.37 * t) / 0.37, 1e-13));
    CHECK(c.c12 * c.c12 <= c.c11 * c.c22);
    CHECK(kModel.covariance(0.0).c11 == 0.0);

    // Scale s on an interval multiplies that interval's contribution by s².
    const G2pp scaled = kModel.withScale({1.0, 2.0}, {1.0, 2.0});
    const G2pp::Covariance a = kModel.covariance(1.0), b = kModel.covariance(2.0), d = kModel.covariance(3.0);
    CHECK_THAT(scaled.covariance(1.0).c11, WithinRel(a.c11, 1e-13));
    CHECK_THAT(scaled.covariance(2.0).c12, WithinRel(a.c12 + 4.0 * (b.c12 - a.c12), 1e-13));
    CHECK_THAT(scaled.covariance(3.0).c22, WithinRel(a.c22 + 4.0 * (d.c22 - a.c22), 1e-13)); // the last scale continues
    CHECK_THAT(kModel.withScale({5.0}, {1.0}).covariance(7.0).c11, WithinRel(kModel.covariance(7.0).c11, 1e-13));

    CHECK_THROWS_AS(G2pp(0.1, 0.1, 0.0, 0.01, 0.0), rl::InputError);
    CHECK_THROWS_AS(G2pp(0.1, 0.1, 0.01, 0.01, 1.5), rl::InputError);
    CHECK_THROWS_AS(G2pp(0.1, 0.1, 0.01, 0.01, 0.0, {2.0, 1.0}, {1.0, 1.0}), rl::InputError);
}

TEST_CASE("G2++ fits the curve: the deflated bond is a martingale for every maturity", "[models][g2pp]") {
    const rl::math::GaussHermite gh(40);
    for (double t : {0.5, 2.0, 10.0}) {
        const G2pp::Covariance c = kModel.covariance(t);
        const double sd1 = std::sqrt(c.c11), beta = c.c12 / c.c11, sd2 = std::sqrt(c.c22 - c.c12 * beta);
        for (double T : {t, t + 1.0, t + 20.0}) {
            const double df = std::exp(-0.035 * T);
            const double mean = gh.integrate([&](double u) { return gh.integrate([&](double v) { return kModel.deflatedBond(df, T, c, sd1 * u, beta * sd1 * u + sd2 * v); }); });
            CHECK_THAT(mean, WithinRel(df, 1e-12));
        }
    }
}

TEST_CASE("G2++ zero-bond option reproduces Brigo-Mercurio's formula", "[models][g2pp]") {
    // Σ² = σ²/(2a³)(1−e^{−a(S−T)})²(1−e^{−2aT}) + η²/(2b³)(1−e^{−b(S−T)})²(1−e^{−2bT})
    //      + 2ρση/(ab(a+b)) (1−e^{−a(S−T)})(1−e^{−b(S−T)})(1−e^{−(a+b)T})          (their eq. 4.10)
    const double a = 0.35, b = 0.02, sigma = 0.011, eta = 0.009, rho = -0.75, T = 2.0, S = 7.0, k = 0.8;
    const double ea = 1.0 - std::exp(-a * (S - T)), eb = 1.0 - std::exp(-b * (S - T));
    const double var = sigma * sigma / (2 * a * a * a) * ea * ea * (1.0 - std::exp(-2 * a * T)) + eta * eta / (2 * b * b * b) * eb * eb * (1.0 - std::exp(-2 * b * T)) +
                       2 * rho * sigma * eta / (a * b * (a + b)) * ea * eb * (1.0 - std::exp(-(a + b) * T));
    const double v = std::sqrt(var), dT = std::exp(-0.03 * T), dS = std::exp(-0.034 * S);
    const double d1 = std::log(dT / (k * dS)) / v + 0.5 * v;
    const double bm = dT * rl::math::normalCdf(d1) - k * dS * rl::math::normalCdf(d1 - v);
    CHECK_THAT(g2ppZeroBondPut(kModel, T, S, dT, dS, k), WithinRel(bm, 1e-12));
    // The same option as a two-bond payoff through the swaption pricer.
    const std::vector<BondFlow> flows{{T, 1.0, dT}, {S, -k, dS}};
    CHECK_THAT(g2ppOptionValue(kModel, T, flows), WithinRel(bm, 1e-10));
}

TEST_CASE("G2++ swaption: semi-closed form equals two-dimensional integration of the payoff", "[models][g2pp]") {
    for (const G2pp& model : {kModel, G2pp(0.05, 0.6, 0.008, 0.012, 0.4), kModel.withScale({2.0, 5.0}, {1.3, 0.7}), G2pp(0.8, 0.01, 0.02, 0.012, -0.97),
                              G2pp(0.1, 0.12, 0.006, 0.005, 0.98)}) {
        for (double strike : {0.02, 0.035, 0.05}) {
            for (bool payer : {true, false}) {
                const auto rights = rights10nc1(strike, payer);
                for (std::size_t i : {std::size_t(0), std::size_t(4), std::size_t(8)}) {
                    const double closed = g2ppSwaptionPrice(model, rights[i]);
                    CHECK_THAT(closed, WithinAbs(g2ppOptionValueByIntegration(model, rights[i].time, rights[i].flows, 1201), 2e-7));
                    CHECK(closed >= std::max(rights[i].forwardValue(), 0.0) - 1e-12);
                }
            }
        }
    }
}

TEST_CASE("G2++ swaption: quadrature has converged; put-call parity holds", "[models][g2pp]") {
    const auto pay = rights10nc1(0.033, true), rec = rights10nc1(0.033, false);
    for (std::size_t i = 0; i < pay.size(); ++i) {
        const double p = g2ppSwaptionPrice(kModel, pay[i]);
        CHECK_THAT(g2ppOptionValue(kModel, pay[i].time, pay[i].flows, 24), WithinAbs(p, 1e-11));
        CHECK_THAT(g2ppOptionValue(kModel, pay[i].time, pay[i].flows, 96), WithinAbs(p, 1e-12));
        CHECK_THAT(p - g2ppSwaptionPrice(kModel, rec[i]), WithinAbs(pay[i].forwardValue(), 1e-11));
    }
}

TEST_CASE("G2++ collapses to the one-factor model when the second factor adds nothing", "[models][g2pp]") {
    const double kappa = 0.04, sigma = 0.009;
    const Lgm1F one = Lgm1F::constantSigma(kappa, sigma, 40.0);
    for (const ExerciseRight& r : rights10nc1(0.035)) {
        const double ref = lgmSwaptionPrice(one, r);
        // No second factor.
        CHECK_THAT(g2ppSwaptionPrice(G2pp(kappa, 0.5, sigma, 0.0, 0.3), r), WithinAbs(ref, 1e-13));
        // Two perfectly correlated factors of the same speed are one factor with the summed volatility.
        CHECK_THAT(g2ppSwaptionPrice(G2pp(kappa, kappa, 0.6 * sigma, 0.4 * sigma, 1.0), r), WithinAbs(ref, 1e-13));
        // A nearly perfectly correlated pair goes through the two-dimensional branch and lands on the same price.
        CHECK_THAT(g2ppSwaptionPrice(G2pp(kappa, kappa, 0.6 * sigma, 0.4 * sigma, 1.0 - 1e-9), r), WithinAbs(ref, 1e-9));
    }
}

TEST_CASE("G2++ decorrelates rates: forward-rate correlation is 1 on the diagonal and falls with distance", "[models][g2pp]") {
    CHECK_THAT(kModel.forwardRateCorrelation(5.0, 5.0), WithinAbs(1.0, 1e-14));
    CHECK(kModel.forwardRateCorrelation(1.0, 2.0) > kModel.forwardRateCorrelation(1.0, 10.0));
    CHECK(kModel.forwardRateCorrelation(1.0, 10.0) < 0.95);
    CHECK_THAT(kModel.forwardRateCorrelation(2.0, 9.0), WithinAbs(kModel.forwardRateCorrelation(9.0, 2.0), 1e-15));
    // One factor: perfectly correlated whatever the maturities.
    const G2pp one(0.04, 0.5, 0.009, 0.0, 0.0);
    CHECK_THAT(one.forwardRateCorrelation(0.5, 30.0), WithinAbs(1.0, 1e-14));
    // A humped volatility term structure needs the two factors to offset at the short end.
    const G2pp humped(0.5, 0.02, 0.012, 0.010, -0.9);
    CHECK(humped.forwardRateVolatility(3.0) > humped.forwardRateVolatility(0.0));
    CHECK(humped.forwardRateVolatility(3.0) > humped.forwardRateVolatility(30.0));
}
