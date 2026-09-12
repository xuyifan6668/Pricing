#include "rl/math/rootfind.hpp"

#include "rl/core/error.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>

using namespace rl::math;
using Catch::Matchers::WithinAbs;

TEST_CASE("Brent: sqrt(2), transcendental, endpoint root", "[math][rootfind]") {
    // The default xTol is 1e-12 (absolute), which is what a rate solve needs; ask for 1e-15 to check
    // that the algorithm itself reaches machine precision.
    RootOptions tight;
    tight.xTol = 1e-15;
    auto r = brent([](double x) { return x * x - 2.0; }, 0.0, 2.0, tight);
    CHECK(r.converged);
    CHECK_THAT(r.x, WithinAbs(std::sqrt(2.0), 1e-14));
    CHECK(r.iterations < 20);

    auto s = brent([](double x) { return std::exp(x) - x - 2.0; }, 0.0, 3.0, tight);
    CHECK_THAT(s.x, WithinAbs(1.1461932206205825, 1e-13));

    auto e = brent([](double x) { return x - 1.0; }, 1.0, 5.0); // root at the left endpoint
    CHECK_THAT(e.x, WithinAbs(1.0, 1e-14));

    auto flat = brent([](double x) { return std::tanh(50 * (x - 0.3)); }, -1.0, 1.0); // steep, tests bisection fallback
    CHECK_THAT(flat.x, WithinAbs(0.3, 1e-11));
}

TEST_CASE("Brent: tolerances and failure modes", "[math][rootfind]") {
    RootOptions opt;
    opt.xTol = 1e-6;
    auto r = brent([](double x) { return x * x - 2.0; }, 0.0, 2.0, opt);
    CHECK_THAT(r.x, WithinAbs(std::sqrt(2.0), 1e-6));
    CHECK_THROWS_AS(brent([](double x) { return x * x + 1.0; }, 0.0, 2.0), rl::InputError); // not bracketed
    opt.maxIter = 2;
    opt.xTol = 1e-15;
    CHECK_THROWS_AS(brent([](double x) { return std::tanh(50 * (x - 0.3)); }, -1.0, 1.0, opt), rl::NumericalError);
    // fTol stops early.
    RootOptions ft;
    ft.fTol = 1e-3;
    auto q = brent([](double x) { return x * x - 2.0; }, 0.0, 2.0, ft);
    CHECK(std::fabs(q.f) <= 1e-3);
}

TEST_CASE("Safeguarded Newton", "[math][rootfind]") {
    auto fdf = [](double x) { return std::make_pair(x * x - 2.0, 2 * x); };
    auto r = newtonSafeguarded(fdf, 1.0, 0.0, 2.0);
    CHECK(r.converged);
    CHECK_THAT(r.x, WithinAbs(std::sqrt(2.0), 1e-14));
    CHECK(r.iterations < 10);
    // Start at the far end of the bracket; derivative zero at the start (x0 = 0) must not blow up.
    auto z = newtonSafeguarded(fdf, 0.0, 0.0, 2.0);
    CHECK_THAT(z.x, WithinAbs(std::sqrt(2.0), 1e-14));
    auto h = newtonSafeguarded(fdf, 2.0, 0.0, 2.0);
    CHECK_THAT(h.x, WithinAbs(std::sqrt(2.0), 1e-14));
    // A function where pure Newton diverges from x0 (f = arctan): safeguard keeps it in the bracket.
    auto at = [](double x) { return std::make_pair(std::atan(x - 0.5), 1.0 / (1.0 + (x - 0.5) * (x - 0.5))); };
    auto a = newtonSafeguarded(at, 20.0, -30.0, 30.0);
    CHECK_THAT(a.x, WithinAbs(0.5, 1e-13));
    CHECK_THROWS_AS(newtonSafeguarded(fdf, 1.0, 2.0, 3.0), rl::InputError); // no sign change
}

TEST_CASE("expandBracket", "[math][rootfind]") {
    double a = 0.0, b = 1.0;
    CHECK(expandBracket([](double x) { return x - 100.0; }, a, b));
    CHECK(a <= 100.0);
    CHECK(b >= 100.0);
    double c = 0.0, d = 1.0;
    CHECK_FALSE(expandBracket([](double x) { return x * x + 1.0; }, c, d, 10));
}
