#include "rl/math/normal.hpp"

#include "rl/aad/real.hpp"
#include "rl/core/error.hpp"
#include "rl/math/quadrature.hpp"

#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>

using rl::aad::Real;
using namespace rl::math;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

TEMPLATE_TEST_CASE("pdf: value at 0, symmetry, integrates to one", "[math][normal]", double, Real) {
    using T = TestType;
    CHECK(value(normalPdf(T(0.0))) == 0.3989422804014327);
    for (double x : {0.3, 1.0, 2.5, 6.0}) CHECK(value(normalPdf(T(x))) == value(normalPdf(T(-x))));
    CHECK_THAT(value(normalPdf(T(1.0))), WithinRel(0.24197072451914337, 1e-15));
    CHECK(value(normalPdf(T(40.0))) == 0.0); // underflow is fine
}

TEMPLATE_TEST_CASE("cdf: goldens, symmetry, tails", "[math][normal]", double, Real) {
    using T = TestType;
    CHECK(value(normalCdf(T(0.0))) == 0.5);
    CHECK_THAT(value(normalCdf(T(-1.959963984540054))), WithinAbs(0.025, 1e-15));
    CHECK_THAT(value(normalCdf(T(1.0))), WithinRel(0.8413447460685429, 1e-15));
    CHECK_THAT(value(normalCdf(T(-1.0))), WithinRel(0.15865525393145707, 1e-15));
    for (double x : {0.1, 0.7, 1.3, 2.2, 3.9, 5.5}) CHECK_THAT(value(normalCdf(T(x))) + value(normalCdf(T(-x))), WithinAbs(1.0, 1e-16));
    // Deep lower tail keeps relative accuracy (this is the erfc-vs-erf point).
    CHECK_THAT(value(normalCdf(T(-8.0))), WithinRel(6.220960574271786e-16, 1e-13));
    CHECK_THAT(value(normalCdf(T(-20.0))), WithinRel(2.7536241186062337e-89, 1e-12));
    CHECK(value(normalCdf(T(40.0))) == 1.0);
    CHECK(value(normalCdf(T(-40.0))) == 0.0);
}

TEMPLATE_TEST_CASE("inverse cdf: goldens and antisymmetry", "[math][normal]", double, Real) {
    using T = TestType;
    CHECK(value(normalInverseCdf(T(0.5))) == 0.0);
    CHECK_THAT(value(normalInverseCdf(T(0.975))), WithinAbs(1.959963984540054, 2e-15));
    CHECK_THAT(value(normalInverseCdf(T(0.8413447460685429))), WithinAbs(1.0, 2e-15));
    CHECK_THAT(value(normalInverseCdf(T(1e-10))), WithinRel(-6.361340902404056, 1e-13));
    for (double p : {0.01, 0.1, 0.3, 0.45}) CHECK_THAT(value(normalInverseCdf(T(p))), WithinAbs(-value(normalInverseCdf(T(1.0 - p))), 1e-14));
    CHECK_THROWS_AS(normalInverseCdf(T(0.0)), rl::InputError);
    CHECK_THROWS_AS(normalInverseCdf(T(1.0)), rl::InputError);
    CHECK_THROWS_AS(normalInverseCdf(T(-0.1)), rl::InputError);
}

TEST_CASE("AS241 round trip cdf(inv(p)) = p to the conditioning limit over [1e-300, 1 - 1e-16]", "[math][normal]") {
    // Log-spaced p in the lower tail; the upper tail is checked by symmetry and near 1 explicitly.
    // Conditioning: one ulp of x moves Φ(x) by |x|·φ(x)/Φ(x)·ulp ≈ x²·ε relative in the tail, so the
    // attainable relative accuracy is ~ε(1 + x²), not ε. That is the tolerance, with a 4× margin for erfc.
    for (double e = -300.0; e <= -0.31; e += 0.37) {
        double p = std::pow(10.0, e);
        double x = normalInverseCdf<double>(p);
        CHECK_THAT(normalCdf<double>(x), WithinRel(p, 1e-15 * (1.0 + x * x)));
    }
    for (double p = 0.001; p < 0.9995; p += 0.00173) {
        double x = normalInverseCdf<double>(p);
        CHECK_THAT(normalCdf<double>(x), WithinAbs(p, 1e-15));
    }
    // Near one: p = 1 − q is only representable to 1.1e-16 absolute, so the complement can only be
    // recovered to that absolute accuracy (the upper tail is where the symmetry test carries the load).
    for (double q : {1e-12, 1e-8, 1e-4}) {
        double p = 1.0 - q;
        double x = normalInverseCdf<double>(p);
        CHECK_THAT(1.0 - normalCdf<double>(x), WithinAbs(q, 4e-16));
    }
    // Far tail golden (mpmath, 30 digits: −37.04709629936119923722).
    CHECK_THAT(normalInverseCdf<double>(1e-300), WithinRel(-37.0470962993612, 1e-13));
}

TEST_CASE("pdf integrates to one and is the derivative of the cdf", "[math][normal]") {
    GaussHermite gh(40);
    // E[ pdf(Z) / pdf(Z) ] trivially 1; instead check ∫ pdf = 1 by Simpson on [-10, 10].
    const int n = 4000;
    const double a = -10.0, b = 10.0, h = (b - a) / n;
    double s = normalPdf<double>(a) + normalPdf<double>(b);
    for (int i = 1; i < n; ++i) s += (i % 2 ? 4.0 : 2.0) * normalPdf<double>(a + i * h);
    CHECK_THAT(s * h / 3.0, WithinAbs(1.0, 1e-12));
    for (double x : {-2.0, -0.5, 0.0, 0.8, 3.0}) {
        double d = 1e-5;
        CHECK_THAT((normalCdf<double>(x + d) - normalCdf<double>(x - d)) / (2 * d), WithinAbs(normalPdf<double>(x), 1e-9));
    }
    // Gauss–Hermite: E[Z²] = 1 as a smoke test that the two modules agree on the normalisation.
    CHECK_THAT(gh.integrate([](double z) { return z * z; }), WithinAbs(1.0, 1e-12));
}

TEST_CASE("Real and double agree bit for bit", "[math][normal][aad]") {
    for (double x : {-3.0, -0.2, 0.0, 0.9, 4.5}) {
        Real xr = Real::input(x);
        CHECK(normalPdf(xr).value() == normalPdf<double>(x));
        CHECK(normalCdf(xr).value() == normalCdf<double>(x));
        CHECK(normalCdf(xr).isActive());
    }
    for (double p : {0.001, 0.2, 0.5, 0.77, 0.9999}) {
        Real pr = Real::input(p);
        CHECK(normalInverseCdf(pr).value() == normalInverseCdf<double>(p));
    }
    rl::aad::Tape::current().reset();
}
