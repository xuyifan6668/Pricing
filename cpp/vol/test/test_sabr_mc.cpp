// The SABR model by Monte Carlo: exact where the model is Bachelier, agrees with Hagan where the
// expansion is known to be accurate, and is reproducible.
#include "rl/vol/sabr_mc.hpp"

#include "rl/core/error.hpp"
#include "rl/vol/sabr.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>

using namespace rl;
using namespace rl::vol;
using Catch::Matchers::WithinAbs;

namespace {
const double kBp = 1e-4;
market::SabrParams<double> params(double alpha, double beta, double rho, double nu) {
    market::SabrParams<double> p;
    p.alpha = alpha;
    p.beta = beta;
    p.rho = rho;
    p.nu = nu;
    return p;
}
} // namespace

TEST_CASE("SABR Monte Carlo: with beta 0 and no vol of vol the model is Bachelier", "[vol][sabr][mc]") {
    const double F = 0.035, T = 5.0, sigma = 0.0085;
    const SabrMcResult r = sabrMonteCarlo(params(sigma, 0.0, 0.0, 0.0), F, T, {F - 0.01, F, F + 0.01}, {.paths = 40000, .stepsPerYear = 4, .seed = 1});
    for (const SabrMcPoint& p : r.points) CHECK_THAT(p.normalVol, WithinAbs(sigma, 4.0 * p.normalVolError + 0.05 * kBp));
    CHECK(r.absorbed == 0.0);
}

TEST_CASE("SABR Monte Carlo: agrees with Hagan's expansion where it is accurate (short expiry, moderate vol of vol)", "[vol][sabr][mc]") {
    const double F = 0.036, T = 1.0;
    const market::SabrParams<double> p = params(0.045, 0.5, -0.2, 0.3);
    std::vector<double> strikes;
    for (int off = -100; off <= 100; off += 50) strikes.push_back(F + off * kBp);
    const SabrMcResult r = sabrMonteCarlo(p, F, T, strikes, {.paths = 100000, .stepsPerYear = 200, .seed = 2});
    for (const SabrMcPoint& pt : r.points) {
        INFO("K " << pt.strike << ": MC " << pt.normalVol / kBp << " ± " << pt.normalVolError / kBp << " bp, Hagan " << pt.haganNormalVol / kBp << " bp");
        CHECK_THAT(pt.normalVol, WithinAbs(pt.haganNormalVol, 4.0 * pt.normalVolError + 0.3 * kBp));
    }
}

TEST_CASE("SABR Monte Carlo: the forward is a martingale up to the Euler error; runs are reproducible", "[vol][sabr][mc]") {
    const market::SabrParams<double> p = params(0.045, 0.5, -0.3, 0.5);
    const SabrMcOptions o{.paths = 20000, .stepsPerYear = 100, .seed = 3};
    const SabrMcResult a = sabrMonteCarlo(p, 0.036, 5.0, {0.036}, o), b = sabrMonteCarlo(p, 0.036, 5.0, {0.036}, o);
    CHECK(a.points[0].price == b.points[0].price);
    CHECK(std::fabs(a.forwardError) < 2.0 * kBp);
    CHECK(a.points[0].price > 0.0);
    CHECK_THROWS_AS(sabrMonteCarlo(p, 0.036, 0.0, {0.036}), rl::InputError);
    CHECK_THROWS_AS(sabrMonteCarlo(p, 0.036, 1.0, {0.036}, {.paths = 1001}), rl::InputError);
}
