#include "rl/math/scalar.hpp"

#include "rl/aad/real.hpp"

#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>

using rl::aad::Real;

static_assert(rl::math::Scalar<double>);
static_assert(rl::math::Scalar<Real>);
static_assert(!rl::math::Scalar<const char*>);

TEMPLATE_TEST_CASE("max/min/positivePart/sign under both scalars", "[math][scalar]", double, Real) {
    using T = TestType;
    T a = T(1.0), b = T(2.0), n = T(-3.0);
    CHECK(rl::math::value(rl::math::max(a, b)) == 2.0);
    CHECK(rl::math::value(rl::math::min(a, b)) == 1.0);
    CHECK(rl::math::value(rl::math::positivePart(n)) == 0.0);
    CHECK(rl::math::value(rl::math::positivePart(b)) == 2.0);
    CHECK(rl::math::sign(n) == -1.0);
    CHECK(rl::math::sign(T(0.0)) == 0.0);
    CHECK(rl::math::isFinite(a));
    volatile double zero = 0.0;
    CHECK_FALSE(rl::math::isFinite(T(1.0) / T(zero)));
}

TEMPLATE_TEST_CASE("ADL pattern: unqualified calls resolve for both scalars", "[math][scalar]", double, Real) {
    using T = TestType;
    using rl::math::exp;
    using rl::math::log;
    using rl::math::sqrt;
    T x = T(0.5);
    T y = exp(log(sqrt(x * x)));
    CHECK(rl::math::value(y) == 0.5);
}
