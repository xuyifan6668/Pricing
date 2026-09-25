// YieldCurve<T>: nodes, interpolation, time convention, extrapolation, the turn-of-year jump layer.
#include "rl/market/yield_curve.hpp"

#include "rl/aad/real.hpp"
#include "rl/core/daycount.hpp"
#include "rl/core/error.hpp"

#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <vector>

using namespace rl::market;
using rl::aad::Real;
using rl::core::Date;
using rl::math::value;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace {
const Date kRef(2026, 9, 4);
const std::vector<double> kTimes{0.25, 0.5, 1.0, 2.0, 5.0, 10.0, 30.0};
const std::vector<double> kZeros{0.038, 0.037, 0.035, 0.033, 0.033, 0.0345, 0.035};

template <class T>
std::vector<T> logDfNodes() { // y_i = r_i t_i
    std::vector<T> ys;
    for (std::size_t i = 0; i < kTimes.size(); ++i) ys.push_back(T(kZeros[i] * kTimes[i]));
    return ys;
}
} // namespace

TEMPLATE_TEST_CASE("YieldCurve: df(0) = 1, nodes reproduced, positivity, both interpolations", "[market][curve]", double, Real) {
    using T = TestType;
    for (auto kind : {CurveInterpolation::LinearLogDf, CurveInterpolation::MonotoneConvex}) {
        YieldCurve<T> c(kRef, "USD-SOFR-OIS", kTimes, logDfNodes<T>(), kind);
        CHECK(c.referenceDate() == kRef);
        CHECK(c.name() == "USD-SOFR-OIS");
        CHECK(c.interpolation() == kind);
        CHECK(value(c.df(0.0)) == 1.0);
        CHECK(value(c.df(kRef)) == 1.0);
        CHECK(c.nodeTimes() == kTimes);
        REQUIRE(c.nodeValues().size() == kTimes.size());
        for (std::size_t i = 0; i < kTimes.size(); ++i) {
            CHECK_THAT(value(c.df(kTimes[i])), WithinRel(std::exp(-kZeros[i] * kTimes[i]), 1e-13));
            CHECK(value(c.nodeValues()[i]) == kZeros[i] * kTimes[i]);
        }
        for (double t = 0.0; t <= 40.0; t += 0.173) CHECK(value(c.df(t)) > 0.0);
    }
}

TEST_CASE("YieldCurve: time() is ACT/365F from the reference date", "[market][curve]") {
    YieldCurve<double> c(kRef, "c", kTimes, logDfNodes<double>(), CurveInterpolation::LinearLogDf);
    CHECK(c.time(kRef) == 0.0);
    CHECK_THAT(c.time(Date(2027, 9, 4)), WithinAbs(365.0 / 365.0, 1e-15));
    CHECK_THAT(c.time(Date(2026, 9, 5)), WithinAbs(1.0 / 365.0, 1e-15));
    CHECK(c.time(Date(2026, 9, 3)) < 0.0);
    CHECK(c.df(Date(2027, 9, 4)) == c.df(c.time(Date(2027, 9, 4))));
}

TEST_CASE("YieldCurve: linear on log-df gives piecewise-flat forwards and flat-forward extrapolation", "[market][curve]") {
    YieldCurve<double> c(kRef, "c", kTimes, logDfNodes<double>(), CurveInterpolation::LinearLogDf);
    // Between nodes the continuously-compounded forward is constant: log(P(t1)/P(t2))/(t2−t1).
    auto fwd = [&](double a, double b) { return std::log(c.df(a) / c.df(b)) / (b - a); };
    CHECK_THAT(fwd(1.0, 1.5), WithinRel(fwd(1.5, 2.0), 1e-12));
    CHECK_THAT(fwd(1.2, 1.3), WithinRel(fwd(1.0, 2.0), 1e-12));
    // Before the first node: flat at the first interval's forward (y(0) = 0 anchor).
    CHECK_THAT(fwd(0.0, 0.1), WithinRel(fwd(0.0, 0.25), 1e-12));
    // Beyond the last node: flat at the last interval's forward.
    const double fLast = fwd(10.0, 30.0);
    CHECK_THAT(fwd(30.0, 35.0), WithinRel(fLast, 1e-12));
    CHECK_THAT(fwd(40.0, 41.0), WithinRel(fLast, 1e-12));
}

TEST_CASE("YieldCurve: withNodeValues keeps everything but the ordinates", "[market][curve]") {
    YieldCurve<double> c(kRef, "c", kTimes, logDfNodes<double>(), CurveInterpolation::MonotoneConvex);
    std::vector<double> ys = c.nodeValues();
    ys[3] += 0.001;
    YieldCurve<double> d = c.withNodeValues(ys);
    CHECK(d.referenceDate() == c.referenceDate());
    CHECK(d.name() == c.name());
    CHECK(d.interpolation() == c.interpolation());
    CHECK(d.nodeTimes() == c.nodeTimes());
    CHECK(d.nodeValues()[3] == ys[3]);
    CHECK(c.nodeValues()[3] != ys[3]); // original untouched
    CHECK(d.df(2.0) < c.df(2.0));
    CHECK(d.df(0.25) == c.df(0.25)); // monotone convex is local: the 3m node does not see a 2y bump
    CHECK_THROWS_AS(c.withNodeValues(std::vector<double>(3, 0.1)), rl::InputError);
    CHECK_THROWS_AS(c.withNodeValues(std::vector<double>(kTimes.size(), std::nan(""))), rl::InputError);
}

TEST_CASE("YieldCurve: turn-of-year jump layer moves the forward over the window and nothing else", "[market][curve][turn]") {
    using C = YieldCurve<double>;
    const double a = 0.32, b = 0.33; // ~4 days
    C plain(kRef, "c", kTimes, logDfNodes<double>(), CurveInterpolation::LinearLogDf);
    C turned(kRef, "c", kTimes, logDfNodes<double>(), CurveInterpolation::LinearLogDf, {TurnJump<double>{a, b, 8e-4}});
    REQUIRE(turned.turns().size() == 1);
    CHECK(turned.df(0.0) == 1.0);
    // Before the window: identical.
    for (double t : {0.05, 0.2, 0.32}) CHECK(turned.df(t) == plain.df(t));
    // Over the window the continuously-compounded forward is 8 bp higher.
    CHECK_THAT(std::log(turned.df(a) / turned.df(b)) / (b - a), WithinAbs(std::log(plain.df(a) / plain.df(b)) / (b - a) + 8e-4, 1e-12));
    // After it: a constant factor exp(−J·(b−a)), so forwards beyond the window are unchanged.
    for (double t : {0.5, 2.0, 10.0}) CHECK_THAT(turned.df(t) / plain.df(t), WithinRel(std::exp(-8e-4 * (b - a)), 1e-13));
    CHECK_THAT(std::log(turned.df(1.0) / turned.df(2.0)), WithinRel(std::log(plain.df(1.0) / plain.df(2.0)), 1e-12));
    CHECK(turned.dfSmooth(5.0) == plain.df(5.0));
    // withNodeValues keeps the turn; withTurnJumps replaces it.
    CHECK(turned.withNodeValues(logDfNodes<double>()).turns().size() == 1);
    CHECK(turned.withTurnJumps({0.0}).df(10.0) == plain.df(10.0));
    CHECK_THROWS_AS(turned.withTurnJumps({0.0, 0.0}), rl::InputError);
    CHECK_THROWS_AS(C(kRef, "c", kTimes, logDfNodes<double>(), CurveInterpolation::LinearLogDf, {TurnJump<double>{0.5, 0.4, 1e-4}}), rl::InputError);
}

TEST_CASE("YieldCurve validation", "[market][curve]") {
    using C = YieldCurve<double>;
    const auto k = CurveInterpolation::LinearLogDf;
    CHECK_THROWS_AS(C(kRef, "c", {}, {}, k), rl::InputError);
    CHECK_THROWS_AS(C(kRef, "c", {0.0, 1.0}, {0.0, 0.03}, k), rl::InputError);   // t = 0 is not a free node
    CHECK_THROWS_AS(C(kRef, "c", {1.0, 1.0}, {0.03, 0.03}, k), rl::InputError);   // not strictly increasing
    CHECK_THROWS_AS(C(kRef, "c", {2.0, 1.0}, {0.03, 0.03}, k), rl::InputError);
    CHECK_THROWS_AS(C(kRef, "c", {1.0, 2.0}, {0.03}, k), rl::InputError);         // size mismatch
    CHECK_THROWS_AS(C(kRef, "c", {1.0, 2.0}, {0.03, std::nan("")}, k), rl::InputError);
    CHECK_NOTHROW(C(kRef, "c", {1.0}, {0.03}, k));                                // one node is a flat curve
    C flat(kRef, "c", {1.0}, {0.03}, k);
    CHECK_THAT(flat.df(2.0), WithinRel(std::exp(-0.06), 1e-14));
}

TEST_CASE("YieldCurve: Real and double agree bit for bit; Real nodes stay active", "[market][curve][aad]") {
    std::vector<Real> yr;
    for (double y : logDfNodes<double>()) yr.push_back(Real::input(y));
    for (auto kind : {CurveInterpolation::LinearLogDf, CurveInterpolation::MonotoneConvex}) {
        YieldCurve<double> cd(kRef, "c", kTimes, logDfNodes<double>(), kind);
        YieldCurve<Real> cr(kRef, "c", kTimes, yr, kind);
        for (double t = 0.05; t <= 35.0; t += 0.37) {
            CHECK(cd.df(t) == cr.df(t).value());
            CHECK(cr.df(t).isActive());
        }
        CHECK_FALSE(cr.df(0.0).isActive()); // exactly 1, a constant
    }
    rl::aad::Tape::current().reset();
}
