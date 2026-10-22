// European swaptions under LGM, tested on what a swaption is rather than on stored numbers.
#include "rl/models/lgm_swaption.hpp"

#include "rl/core/error.hpp"
#include "rl/math/normal.hpp"
#include "test_support.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>

using namespace rl::models;
using namespace rl::models::test;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

TEST_CASE("Jamshidian closed form equals direct integration of the payoff", "[models][swaption]") {
    for (double kappa : {0.0, 0.03, 0.10}) {
        const Lgm1F model = Lgm1F::constantSigma(kappa, 0.009, 30.0);
        for (double strike : {0.02, 0.035, 0.05}) {
            for (bool payer : {true, false}) {
                for (const ExerciseRight& r : rights10nc1(strike, payer)) {
                    const double zeta = model.zeta(r.time);
                    const double closed = lgmOptionValue(model, zeta, r.flows);
                    // Simpson is second order across the payoff's kink at the boundary, so use a fine grid.
                    const double integral = lgmOptionValueByIntegration(model, zeta, r.flows, 40001);
                    CHECK_THAT(closed, WithinAbs(integral, 1e-8));
                    CHECK(closed >= 0.0);
                }
            }
        }
    }
}

TEST_CASE("Put-call parity: payer - receiver = value of the forward swap, for any model parameters", "[models][swaption]") {
    const Lgm1F model = Lgm1F::constantSigma(0.05, 0.012, 30.0);
    const auto pay = rights10nc1(0.033, true), rec = rights10nc1(0.033, false);
    for (std::size_t i = 0; i < pay.size(); ++i)
        CHECK_THAT(lgmSwaptionPrice(model, pay[i]) - lgmSwaptionPrice(model, rec[i]), WithinAbs(pay[i].forwardValue(), 1e-14));
}

TEST_CASE("Limits: no variance gives intrinsic; value rises with variance; never below intrinsic", "[models][swaption]") {
    const Lgm1F shape(0.03, {1.0}, {1.0});
    for (double strike : {0.02, 0.05}) {
        const ExerciseRight r = rights10nc1(strike)[3];
        const double intrinsic = std::max(r.forwardValue(), 0.0);
        CHECK_THAT(lgmOptionValue(shape, 1e-14, r.flows), WithinAbs(intrinsic, 1e-9));
        double prev = intrinsic;
        for (double zeta : {1e-6, 1e-5, 1e-4, 1e-3, 1e-2}) {
            const double v = lgmOptionValue(shape, zeta, r.flows);
            CHECK(v >= prev); // far from the money the time value of a tiny variance underflows
            prev = v;
        }
        CHECK(prev > intrinsic + 1e-4);
    }
}

TEST_CASE("One-period swaption is a zero-bond option: matches the Hull-White bond option formula", "[models][swaption][hw]") {
    // Right to receive 1 at s = t and pay (1 + K tau) at e. Under Hull-White with constant sigma the
    // price is D_s N(d1) - (1 + K tau) D_e N(d2) with
    //   v = sigma/kappa (1 - exp(-kappa (e - t))) sqrt((1 - exp(-2 kappa t))/(2 kappa)).
    const double kappa = 0.04, sigma = 0.01, t = 2.0, e = 3.0, k = 0.035, tau = 1.0;
    const Lgm1F model = Lgm1F::constantSigma(kappa, sigma, 30.0);
    const double ds = std::exp(-0.03 * t), de = std::exp(-0.032 * e);
    const std::vector<BondFlow> flows{{t, 1.0, ds}, {e, -(1.0 + k * tau), de}};
    const double v = sigma / kappa * (1.0 - std::exp(-kappa * (e - t))) * std::sqrt((1.0 - std::exp(-2.0 * kappa * t)) / (2.0 * kappa));
    const double d1 = std::log(ds / ((1.0 + k * tau) * de)) / v + 0.5 * v, d2 = d1 - v;
    const double hw = ds * rl::math::normalCdf(d1) - (1.0 + k * tau) * de * rl::math::normalCdf(d2);
    CHECK_THAT(lgmOptionValue(model, model.zeta(t), flows), WithinRel(hw, 1e-11));
}

TEST_CASE("The exercise boundary: payers exercise when rates are high, receivers when low, at the same state", "[models][swaption]") {
    const Lgm1F model = Lgm1F::constantSigma(0.03, 0.009, 30.0);
    const ExerciseRight pay = rights10nc1(0.035, true)[2], rec = rights10nc1(0.035, false)[2];
    const double zeta = model.zeta(pay.time);
    const ExerciseBoundary bp = exerciseBoundary(model, zeta, pay.flows), br = exerciseBoundary(model, zeta, rec.flows);
    REQUIRE(bp.exists);
    REQUIRE(br.exists);
    CHECK(bp.exerciseAbove);
    CHECK_FALSE(br.exerciseAbove);
    CHECK_THAT(bp.x, WithinAbs(br.x, 1e-12));
    CHECK_THAT(deflatedValue(model, zeta, pay.flows, bp.x), WithinAbs(0.0, 1e-14));
    // A deeply in-the-money right with almost no variance has no boundary: always exercised.
    const ExerciseRight deep = rights10nc1(0.001, true)[2];
    CHECK_FALSE(exerciseBoundary(model, 1e-10, deep.flows).exists);
    // A payoff that changes sign twice is not a swaption.
    const std::vector<BondFlow> butterfly{{1.0, 1.0, 0.97}, {2.0, -2.0, 0.94}, {3.0, 1.0002, 0.91}};
    CHECK_THROWS_AS(exerciseBoundary(Lgm1F::constantSigma(0.0, 0.01, 10.0), 5e-4, butterfly), rl::NumericalError);
}

TEST_CASE("Market formula: implied vol inverts the price; payer and receiver share one vol", "[models][swaption]") {
    const ExerciseRight pay = rights10nc1(0.04, true)[4], rec = rights10nc1(0.04, false)[4];
    const double sigma = 0.0085;
    CHECK_THAT(impliedNormalVol(pay, marketSwaptionPrice(pay, sigma)), WithinRel(sigma, 1e-11));
    CHECK_THAT(impliedNormalVol(rec, marketSwaptionPrice(rec, sigma)), WithinRel(sigma, 1e-11));
    CHECK_THAT(marketSwaptionPrice(pay, sigma) - marketSwaptionPrice(rec, sigma), WithinAbs(pay.forwardValue(), 1e-15));
}
