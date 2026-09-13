#include "rl/math/interpolation.hpp"

#include "rl/aad/real.hpp"
#include "rl/core/error.hpp"

#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <vector>

using rl::aad::Real;
using namespace rl::math;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace {
template <class T>
std::vector<T> lift(const std::vector<double>& v) {
    return std::vector<T>(v.begin(), v.end());
}
// Simpson on [a, b] with 2000 panels: enough for the C¹ interpolants here.
template <class F>
double simpson(F&& f, double a, double b) {
    const int n = 2000;
    const double h = (b - a) / n;
    double s = f(a) + f(b);
    for (int i = 1; i < n; ++i) s += (i % 2 ? 4.0 : 2.0) * f(a + i * h);
    return s * h / 3.0;
}
} // namespace

TEMPLATE_TEST_CASE("Linear: nodes, midpoints, slopes, extrapolation", "[math][interp]", double, Real) {
    using T = TestType;
    std::vector<double> xs{0.0, 1.0, 2.0, 4.0};
    std::vector<double> ys{0.0, 0.03, 0.05, 0.11};
    Linear<T> lin(xs, lift<T>(ys));
    for (std::size_t i = 0; i < xs.size(); ++i) CHECK(value(lin(xs[i])) == ys[i]);
    CHECK_THAT(value(lin(0.5)), WithinAbs(0.015, 1e-16));
    CHECK_THAT(value(lin(3.0)), WithinAbs(0.08, 1e-16));
    CHECK_THAT(value(lin.derivative(0.5)), WithinAbs(0.03, 1e-16));
    CHECK_THAT(value(lin.derivative(3.0)), WithinAbs(0.03, 1e-16));
    CHECK_THAT(value(lin.derivative(1.0)), WithinAbs(0.02, 1e-16)); // right-hand slope at a node
    // Linear extrapolation continues the end segments.
    CHECK_THAT(value(lin(5.0)), WithinAbs(0.14, 1e-15));
    CHECK_THAT(value(lin(-1.0)), WithinAbs(-0.03, 1e-15));
    Linear<T> flat(xs, lift<T>(ys), Extrapolation::Flat);
    CHECK(value(flat(5.0)) == 0.11);
    CHECK(value(flat(-1.0)) == 0.0);
    CHECK(value(flat.derivative(5.0)) == 0.0);
}

TEST_CASE("Linear validation", "[math][interp]") {
    CHECK_THROWS_AS(Linear<double>({0.0}, {1.0}), rl::InputError);
    CHECK_THROWS_AS(Linear<double>({0.0, 1.0}, {1.0}), rl::InputError);
    CHECK_THROWS_AS(Linear<double>({1.0, 1.0}, {1.0, 2.0}), rl::InputError);
    CHECK_THROWS_AS(Linear<double>({1.0, 0.5}, {1.0, 2.0}), rl::InputError);
}

TEMPLATE_TEST_CASE("LogLinearDf: flat forward between nodes", "[math][interp]", double, Real) {
    using T = TestType;
    std::vector<double> ts{0.0, 1.0, 2.0, 5.0};
    std::vector<double> dfs{1.0, std::exp(-0.03), std::exp(-0.03 - 0.035), std::exp(-0.03 - 0.035 - 3 * 0.04)};
    LogLinearDf<T> p(ts, lift<T>(dfs));
    for (std::size_t i = 0; i < ts.size(); ++i) CHECK_THAT(value(p(ts[i])), WithinRel(dfs[i], 1e-15));
    CHECK_THAT(value(p.forward(0.5)), WithinAbs(0.03, 1e-14));
    CHECK_THAT(value(p.forward(1.5)), WithinAbs(0.035, 1e-14));
    CHECK_THAT(value(p.forward(3.7)), WithinAbs(0.04, 1e-14));
    CHECK_THAT(value(p.forward(7.0)), WithinAbs(0.04, 1e-14)); // flat-forward extrapolation
    CHECK_THAT(value(p(7.0)), WithinRel(dfs[3] * std::exp(-2 * 0.04), 1e-14));
    CHECK_THAT(value(p.derivative(1.5)), WithinRel(-0.035 * value(p(1.5)), 1e-13));
    CHECK_THROWS_AS(LogLinearDf<T>({0.0, 1.0}, {T(1.0), T(-0.5)}), rl::InputError);
}

TEMPLATE_TEST_CASE("MonotoneConvex: flat curve is exactly flat", "[math][interp]", double, Real) {
    using T = TestType;
    const double r = 0.04;
    std::vector<double> ts{0.0, 0.5, 1.0, 2.0, 3.0, 5.0, 10.0};
    std::vector<double> ys;
    for (double t : ts) ys.push_back(r * t);
    MonotoneConvex<T> mc(ts, lift<T>(ys));
    for (const auto& f : mc.discreteForwards()) CHECK_THAT(value(f), WithinAbs(r, 1e-15));
    for (const auto& f : mc.nodeForwards()) CHECK_THAT(value(f), WithinAbs(r, 1e-15));
    for (double t : {0.1, 0.75, 1.5, 2.9, 4.0, 7.5, 12.0}) {
        CHECK_THAT(value(mc.forward(t)), WithinAbs(r, 1e-14));
        CHECK_THAT(value(mc(t)), WithinAbs(r * t, 1e-14));
    }
}

TEMPLATE_TEST_CASE("MonotoneConvex: node reproduction, continuity, integral identity, positivity", "[math][interp]", double, Real) {
    using T = TestType;
    // A steep upward-sloping OIS curve: zero rates 2% at 6m rising to 5.2% at 30y.
    std::vector<double> ts{0.0, 0.5, 1.0, 2.0, 3.0, 5.0, 7.0, 10.0, 15.0, 20.0, 30.0};
    std::vector<double> zr{0.0, 0.020, 0.024, 0.030, 0.034, 0.039, 0.042, 0.045, 0.048, 0.050, 0.052};
    std::vector<double> ys;
    for (std::size_t i = 0; i < ts.size(); ++i) ys.push_back(zr[i] * ts[i]);
    MonotoneConvex<T> mc(ts, lift<T>(ys));

    SECTION("nodes reproduced") {
        for (std::size_t i = 0; i < ts.size(); ++i) CHECK_THAT(value(mc(ts[i])), WithinAbs(ys[i], 1e-13));
    }
    SECTION("discrete forwards are the interval averages") {
        for (std::size_t i = 1; i < ts.size(); ++i)
            CHECK_THAT(value(mc.discreteForwards()[i - 1]), WithinAbs((ys[i] - ys[i - 1]) / (ts[i] - ts[i - 1]), 1e-14));
    }
    SECTION("forward is continuous at the nodes and equals the node forward") {
        for (std::size_t i = 1; i + 1 < ts.size(); ++i) {
            double fl = value(mc.forward(ts[i] - 1e-9));
            double fr = value(mc.forward(ts[i] + 1e-9));
            CHECK_THAT(fl, WithinAbs(fr, 1e-6));
            CHECK_THAT(value(mc.forward(ts[i])), WithinAbs(value(mc.nodeForwards()[i]), 1e-13));
        }
    }
    SECTION("y' = f and the integral identity") {
        for (std::size_t i = 1; i < ts.size(); ++i) {
            double integral = simpson([&](double t) { return value(mc.forward(t)); }, ts[i - 1], ts[i]);
            CHECK_THAT(integral, WithinAbs(ys[i] - ys[i - 1], 1e-9));
        }
        for (double t : {0.25, 1.4, 2.5, 4.0, 6.0, 8.8, 12.0, 17.0, 25.0}) {
            double h = 1e-5;
            double num = (value(mc(t + h)) - value(mc(t - h))) / (2 * h);
            CHECK_THAT(num, WithinAbs(value(mc.forward(t)), 1e-8));
            CHECK(value(mc.derivative(t)) == value(mc.forward(t)));
        }
    }
    SECTION("positivity and boundedness of forwards on positive data") {
        for (double t = 0.0; t <= 30.0; t += 0.05) {
            double f = value(mc.forward(t));
            CHECK(f >= 0.0);
            CHECK(f < 0.2);
        }
    }
    SECTION("flat-forward extrapolation beyond the last node") {
        double fN = value(mc.nodeForwards().back());
        CHECK_THAT(value(mc.forward(35.0)), WithinAbs(fN, 1e-14));
        CHECK_THAT(value(mc(35.0)), WithinAbs(ys.back() + 5.0 * fN, 1e-12));
    }
}

TEST_CASE("MonotoneConvex: the four regimes, hand-picked", "[math][interp]") {
    // Three intervals of unit length; choose y so that the middle interval's g0, g1 land in each regime.
    // With interior node forwards from the central formula we cannot dial g0/g1 directly, so instead
    // verify the defining property in each case: integral over the interval equals the discrete forward,
    // and g is continuous. Cases: convex-up, concave, hump, and the plain quadratic regime.
    auto check = [](std::vector<double> ys) {
        std::vector<double> ts{0.0, 1.0, 2.0, 3.0, 4.0};
        MonotoneConvex<double> mc(ts, ys, false); // no collar: exercise the raw regimes
        for (std::size_t i = 1; i < ts.size(); ++i) {
            double integral = simpson([&](double t) { return mc.forward(t); }, ts[i - 1], ts[i]);
            CHECK_THAT(integral, WithinAbs(ys[i] - ys[i - 1], 1e-9));
            CHECK_THAT(mc(ts[i]), WithinAbs(ys[i], 1e-12));
        }
        for (double t = 0.01; t < 4.0; t += 0.01) {
            double fl = mc.forward(t - 5e-7), fr = mc.forward(t + 5e-7);
            CHECK_THAT(fl, WithinAbs(fr, 1e-4)); // continuity everywhere, including at η
        }
    };
    check({0.0, 0.01, 0.04, 0.09, 0.16});   // convex: f^d = 1, 3, 5, 7 (%)
    check({0.0, 0.07, 0.09, 0.10, 0.105});  // concave
    check({0.0, 0.02, 0.07, 0.09, 0.10});   // hump
    check({0.0, 0.03, 0.05, 0.10, 0.11});   // zig-zag: exercises the sign-change branches
}

TEST_CASE("MonotoneConvex validation", "[math][interp]") {
    CHECK_THROWS_AS(MonotoneConvex<double>({0.5, 1.0}, {0.01, 0.02}), rl::InputError); // must start at t = 0
    CHECK_THROWS_AS(MonotoneConvex<double>({0.0, 1.0}, {0.01, 0.02}), rl::InputError); // y(0) must be 0
}

TEMPLATE_TEST_CASE("NaturalCubic: reproduces a cubic with zero end curvature; line gives zero M", "[math][interp]", double, Real) {
    using T = TestType;
    // p(x) = x³ − 3x  has p'' = 6x, which is zero at x = 0 only; use a symmetric-about-zero cubic on
    // [−a, a] with p''(±a) = 0: p(x) = x (a² − x²/3)... p'' = −2x ≠ 0. Simplest exact case: a line.
    std::vector<double> xs{0.0, 1.0, 2.5, 4.0, 6.0};
    std::vector<double> line;
    for (double x : xs) line.push_back(0.2 + 0.1 * x);
    NaturalCubic<T> spl(xs, lift<T>(line));
    for (const auto& m : spl.secondDerivatives()) CHECK_THAT(value(m), WithinAbs(0.0, 1e-14));
    for (double x : {0.3, 1.7, 3.0, 5.5}) {
        CHECK_THAT(value(spl(x)), WithinAbs(0.2 + 0.1 * x, 1e-14));
        CHECK_THAT(value(spl.derivative(x)), WithinAbs(0.1, 1e-14));
        CHECK_THAT(value(spl.secondDerivative(x)), WithinAbs(0.0, 1e-14));
    }
    // Natural condition and node reproduction on curved data.
    std::vector<double> ys{0.30, 0.25, 0.22, 0.21, 0.205};
    NaturalCubic<T> s2(xs, lift<T>(ys));
    CHECK(value(s2.secondDerivatives().front()) == 0.0);
    CHECK(value(s2.secondDerivatives().back()) == 0.0);
    for (std::size_t i = 0; i < xs.size(); ++i) CHECK_THAT(value(s2(xs[i])), WithinAbs(ys[i], 1e-14));
    // C² at interior nodes.
    for (std::size_t i = 1; i + 1 < xs.size(); ++i) {
        CHECK_THAT(value(s2.derivative(xs[i] - 1e-9)), WithinAbs(value(s2.derivative(xs[i] + 1e-9)), 1e-6));
        CHECK_THAT(value(s2.secondDerivative(xs[i] - 1e-9)), WithinAbs(value(s2.secondDerivative(xs[i] + 1e-9)), 1e-6));
    }
    // Flat extrapolation by default.
    CHECK(value(s2(10.0)) == ys.back());
    CHECK(value(s2(-1.0)) == ys.front());
    CHECK_THROWS_AS(NaturalCubic<T>({0.0, 1.0}, {T(1.0), T(2.0)}), rl::InputError);
}

TEST_CASE("NaturalCubic: derivative consistent with finite differences", "[math][interp]") {
    std::vector<double> xs{0.0, 0.5, 1.0, 2.0, 3.0, 5.0};
    std::vector<double> ys{0.80, 0.62, 0.55, 0.50, 0.49, 0.47};
    NaturalCubic<double> s(xs, ys);
    for (double x : {0.2, 0.7, 1.5, 2.5, 4.0}) {
        double h = 1e-5;
        CHECK_THAT((s(x + h) - s(x - h)) / (2 * h), WithinAbs(s.derivative(x), 1e-8));
        CHECK_THAT((s(x + h) - 2 * s(x) + s(x - h)) / (h * h), WithinAbs(s.secondDerivative(x), 1e-4));
    }
}

TEMPLATE_TEST_CASE("Interpolator variant dispatch", "[math][interp]", double, Real) {
    using T = TestType;
    std::vector<double> xs{0.0, 1.0, 2.0};
    Interpolator<T> v = Linear<T>(xs, lift<T>({0.0, 1.0, 4.0}));
    CHECK(value(interpolate(v, 1.5)) == 2.5);
    CHECK(value(interpolateDerivative(v, 1.5)) == 3.0);
    v = NaturalCubic<T>(xs, lift<T>({0.0, 1.0, 4.0}));
    CHECK_THAT(value(interpolate(v, 1.0)), WithinAbs(1.0, 1e-14));
}

TEST_CASE("Real and double interpolation agree bit for bit", "[math][interp][aad]") {
    std::vector<double> ts{0.0, 0.5, 1.0, 2.0, 5.0, 10.0};
    std::vector<double> ys{0.0, 0.012, 0.026, 0.058, 0.17, 0.39};
    std::vector<Real> yr;
    for (double y : ys) yr.push_back(Real::input(y));
    MonotoneConvex<double> md(ts, ys);
    MonotoneConvex<Real> mr(ts, yr);
    NaturalCubic<double> cd(ts, ys);
    NaturalCubic<Real> cr(ts, yr);
    for (double t = 0.0; t <= 12.0; t += 0.37) {
        CHECK(md(t) == mr(t).value());
        CHECK(md.forward(t) == mr.forward(t).value());
        CHECK(cd(t) == cr(t).value());
        CHECK((mr(t).isActive() || t == 0.0));
    }
    rl::aad::Tape::current().reset();
}
