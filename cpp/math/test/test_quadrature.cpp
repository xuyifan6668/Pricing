#include "rl/math/quadrature.hpp"

#include "rl/aad/real.hpp"
#include "rl/core/error.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>

using rl::aad::Real;
using rl::math::GaussHermite;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

TEST_CASE("Gauss-Hermite moments and normalisation", "[math][quadrature]") {
    for (std::size_t n : {1u, 2u, 3u, 5u, 10u, 20u, 64u}) {
        GaussHermite gh(n);
        REQUIRE(gh.size() == n);
        double sw = 0.0;
        for (double w : gh.weights()) {
            CHECK(w > 0.0);
            sw += w;
        }
        CHECK_THAT(sw, WithinAbs(1.0, 1e-14));
        CHECK_THAT(gh.integrate([](double) { return 1.0; }), WithinAbs(1.0, 1e-14));
        CHECK_THAT(gh.integrate([](double z) { return z; }), WithinAbs(0.0, 1e-14));
        if (n >= 2) CHECK_THAT(gh.integrate([](double z) { return z * z; }), WithinAbs(1.0, 1e-13));
        if (n >= 3) CHECK_THAT(gh.integrate([](double z) { return z * z * z * z; }), WithinAbs(3.0, 1e-12));
        if (n >= 4) CHECK_THAT(gh.integrate([](double z) { return std::pow(z, 6); }), WithinAbs(15.0, 1e-11));
        // Sorted ascending, symmetric.
        for (std::size_t i = 1; i < n; ++i) CHECK(gh.nodes()[i] > gh.nodes()[i - 1]);
        for (std::size_t i = 0; i < n; ++i) CHECK(gh.nodes()[i] == -gh.nodes()[n - 1 - i]);
        if (n % 2 == 1) CHECK(gh.nodes()[n / 2] == 0.0);
    }
    GaussHermite one(1);
    CHECK(one.nodes()[0] == 0.0);
    CHECK(one.weights()[0] == 1.0);
}

TEST_CASE("Exactness for degree 2n-1 and failure at 2n", "[math][quadrature]") {
    // E[Z^{2k}] = (2k-1)!!
    auto dfact = [](int k) {
        double p = 1.0;
        for (int i = 2 * k - 1; i > 1; i -= 2) p *= i;
        return p;
    };
    GaussHermite gh(4); // exact up to degree 7
    CHECK_THAT(gh.integrate([](double z) { return std::pow(z, 6); }), WithinRel(dfact(3), 1e-12));
    // Degree 8 is not exact with 4 points (E[Z^8] = 105).
    CHECK(std::fabs(gh.integrate([](double z) { return std::pow(z, 8); }) - 105.0) > 1.0);
    GaussHermite gh5(5);
    CHECK_THAT(gh5.integrate([](double z) { return std::pow(z, 8); }), WithinRel(105.0, 1e-11));
}

TEST_CASE("Lognormal expectation E[exp(sigma Z)] = exp(sigma^2/2)", "[math][quadrature]") {
    GaussHermite gh(20);
    for (double s : {0.1, 0.5, 1.0}) CHECK_THAT(gh.integrate([s](double z) { return std::exp(s * z); }), WithinRel(std::exp(0.5 * s * s), 1e-12));
    // Shifted/scaled form: E[(mu + sigma Z)^2] = mu^2 + sigma^2.
    CHECK_THAT(gh.integrate(0.3, 0.7, [](double x) { return x * x; }), WithinAbs(0.09 + 0.49, 1e-13));
}

TEST_CASE("Integrand returning Real propagates activity", "[math][quadrature][aad]") {
    GaussHermite gh(12);
    Real sigma = Real::input(0.4);
    Real e = gh.integrate([&](double z) { return exp(sigma * z); });
    CHECK(e.isActive());
    CHECK_THAT(e.value(), WithinRel(std::exp(0.08), 1e-12));
    double ed = gh.integrate([](double z) { return std::exp(0.4 * z); });
    CHECK(e.value() == ed);
    rl::aad::Tape::current().reset();
}

TEST_CASE("Order validation", "[math][quadrature]") {
    CHECK_THROWS_AS(GaussHermite(0), rl::InputError);
    CHECK_THROWS_AS(GaussHermite(513), rl::InputError);
}
