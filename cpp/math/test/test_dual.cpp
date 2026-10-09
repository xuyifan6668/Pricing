// Dual<T, N>: forward-mode tangents against closed forms and central differences, under T = double and Real.
#include "rl/math/dual.hpp"

#include "rl/aad/real.hpp"

#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>

using rl::aad::Real;
using namespace rl::math;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

static_assert(Scalar<Dual<double, 1>>);
static_assert(Scalar<Dual<double, 3>>);
static_assert(Scalar<Dual<Real, 3>>);

namespace {
// Uses every operation Dual defines. Written once, evaluated on double (for finite differences) and on Dual.
template <class S>
S composite(const S& x, const S& y) {
    using std::erf;
    using std::erfc;
    using std::exp;
    using std::fabs;
    using std::log;
    using std::pow;
    using std::sqrt;
    return exp(x) * log(y) / sqrt(x * y) + pow(x, y) + erfc(x - y) + pow(y, 2.5) + fabs(x - y) + erf(0.3 * x) + pow(2.0, x) -
           (1.0 - x) / (2.0 + y) + 3.0 / y - y / 4.0;
}
} // namespace

TEMPLATE_TEST_CASE("Dual: product and quotient rules by hand", "[math][dual]", double, Real) {
    using T = TestType;
    using D = Dual<T, 2>;
    const D x = D::variable(T(0.7), 0);
    const D y = D::variable(T(1.9), 1);
    const D f = x * y + x / y - 3.0 * x + 2.0 / y;
    // ∂f/∂x = y + 1/y − 3,  ∂f/∂y = x − x/y² − 2/y²
    CHECK_THAT(f.value(), WithinAbs(0.7 * 1.9 + 0.7 / 1.9 - 2.1 + 2.0 / 1.9, 1e-15));
    CHECK_THAT(value(f.tangent(0)), WithinAbs(1.9 + 1.0 / 1.9 - 3.0, 1e-15));
    CHECK_THAT(value(f.tangent(1)), WithinAbs(0.7 - 0.7 / (1.9 * 1.9) - 2.0 / (1.9 * 1.9), 1e-15));
    D g = x;
    g += y;
    g *= x;
    g -= y;
    g /= y;
    // g = ((x + y)·x − y)/y
    CHECK_THAT(g.value(), WithinAbs(((0.7 + 1.9) * 0.7 - 1.9) / 1.9, 1e-15));
    CHECK_THAT(value(g.tangent(0)), WithinAbs((2.0 * 0.7 + 1.9) / 1.9, 1e-15));
    rl::aad::Tape::current().reset();
}

TEMPLATE_TEST_CASE("Dual: every elementary function matches central differences", "[math][dual]", double, Real) {
    using T = TestType;
    using D = Dual<T, 2>;
    for (auto [x0, y0] : {std::pair{0.7, 1.9}, std::pair{1.3, 0.4}, std::pair{2.2, 2.9}}) {
        const D f = composite(D::variable(T(x0), 0), D::variable(T(y0), 1));
        CHECK_THAT(f.value(), WithinRel(composite(x0, y0), 1e-15));
        const double h = 1e-6;
        const double fx = (composite(x0 + h, y0) - composite(x0 - h, y0)) / (2.0 * h);
        const double fy = (composite(x0, y0 + h) - composite(x0, y0 - h)) / (2.0 * h);
        // Central differences with h = 1e-6 carry ~1e-10 absolute rounding error; 1e-7 relative.
        CHECK_THAT(value(f.tangent(0)), WithinRel(fx, 1e-7));
        CHECK_THAT(value(f.tangent(1)), WithinRel(fy, 1e-7));
    }
    rl::aad::Tape::current().reset();
}

TEST_CASE("Dual: constants carry no tangent; fabs at zero has zero tangent", "[math][dual]") {
    using D = Dual<double, 3>;
    const D c = 2.5;
    for (std::size_t i = 0; i < D::kDirections; ++i) CHECK(c.tangent(i) == 0.0);
    const D x = D::variable(0.0, 1);
    const D a = fabs(x);
    CHECK(a.value() == 0.0);
    CHECK(a.tangent(1) == 0.0);
    CHECK(D::constant(4.0).value() == 4.0);
    // A variable in direction 2 has tangent e_2.
    const D v = D::variable(1.5, 2);
    CHECK(v.tangent(0) == 0.0);
    CHECK(v.tangent(1) == 0.0);
    CHECK(v.tangent(2) == 1.0);
}

TEST_CASE("Dual over Real: bit-identical to the double instantiation, and still active on the tape", "[math][dual][aad]") {
    using DD = Dual<double, 2>;
    using DR = Dual<Real, 2>;
    const DD fd = composite(DD::variable(1.1, 0), DD::variable(0.8, 1));
    const DR fr = composite(DR::variable(Real::input(1.1), 0), DR::variable(Real::input(0.8), 1));
    CHECK(fr.value() == fd.value());
    CHECK(fr.tangent(0).value() == fd.tangent(0));
    CHECK(fr.tangent(1).value() == fd.tangent(1));
    CHECK(fr.base().isActive());
    rl::aad::Tape::current().reset();
}
