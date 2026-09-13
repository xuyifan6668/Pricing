#include "rl/math/tridiagonal.hpp"

#include "rl/aad/real.hpp"
#include "rl/core/error.hpp"
#include "rl/math/rng.hpp"

#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <Eigen/Dense>

using rl::aad::Real;
using rl::math::Tridiagonal;
using Catch::Matchers::WithinAbs;

namespace {
template <class T>
std::vector<T> lift(const std::vector<double>& v) {
    return std::vector<T>(v.begin(), v.end());
}
} // namespace

TEMPLATE_TEST_CASE("Thomas algorithm on a small hand-checkable system", "[math][tridiag]", double, Real) {
    using T = TestType;
    // [ 2 -1  0 ] x = [1]      solution: x = (1, 1, 1)
    // [-1  2 -1 ]     [0]
    // [ 0 -1  2 ]     [1]
    Tridiagonal<T> A(lift<T>({-1.0, -1.0}), lift<T>({2.0, 2.0, 2.0}), lift<T>({-1.0, -1.0}));
    CHECK(A.isDiagonallyDominant());
    auto x = A.solve(lift<T>({1.0, 0.0, 1.0}));
    REQUIRE(x.size() == 3);
    for (const auto& xi : x) CHECK_THAT(rl::math::value(xi), WithinAbs(1.0, 1e-15));
    auto y = A.multiply(x);
    CHECK_THAT(rl::math::value(y[0]), WithinAbs(1.0, 1e-15));
    CHECK_THAT(rl::math::value(y[1]), WithinAbs(0.0, 1e-15));
}

TEST_CASE("Random diagonally dominant systems match Eigen's dense solve", "[math][tridiag]") {
    rl::math::Rng rng(42);
    for (int trial = 0; trial < 20; ++trial) {
        const int n = 5 + trial * 7;
        std::vector<double> a(n - 1), b(n), c(n - 1), d(n);
        for (int i = 0; i < n - 1; ++i) {
            a[i] = rng.uniform01() - 0.5;
            c[i] = rng.uniform01() - 0.5;
        }
        for (int i = 0; i < n; ++i) {
            b[i] = 2.0 + rng.uniform01(); // > |a| + |c| ≤ 1
            d[i] = rng.uniform01() * 10 - 5;
        }
        Tridiagonal<double> A(a, b, c);
        REQUIRE(A.isDiagonallyDominant());
        auto x = A.solve(d);

        Eigen::MatrixXd M = Eigen::MatrixXd::Zero(n, n);
        Eigen::VectorXd rhs(n);
        for (int i = 0; i < n; ++i) {
            M(i, i) = b[i];
            if (i > 0) M(i, i - 1) = a[i - 1];
            if (i + 1 < n) M(i, i + 1) = c[i];
            rhs(i) = d[i];
        }
        Eigen::VectorXd ref = M.colPivHouseholderQr().solve(rhs);
        for (int i = 0; i < n; ++i) CHECK_THAT(x[static_cast<std::size_t>(i)], WithinAbs(ref(i), 1e-12));
        auto back = A.multiply(x);
        for (int i = 0; i < n; ++i) CHECK_THAT(back[static_cast<std::size_t>(i)], WithinAbs(d[static_cast<std::size_t>(i)], 1e-12));
    }
}

TEST_CASE("Diagonal dominance detection and singular pivot", "[math][tridiag]") {
    Tridiagonal<double> weak({1.0, 1.0}, {1.0, 2.0, 1.0}, {1.0, 1.0});
    CHECK_FALSE(weak.isDiagonallyDominant());
    Tridiagonal<double> singular({1.0}, {1.0, 1.0}, {1.0}); // [[1,1],[1,1]]
    CHECK_THROWS_AS(singular.solve({1.0, 2.0}), rl::NumericalError);
    Tridiagonal<double> one({}, {4.0}, {});
    CHECK(one.solve({2.0})[0] == 0.5);
}

TEST_CASE("Validation", "[math][tridiag]") {
    CHECK_THROWS_AS(Tridiagonal<double>({1.0, 1.0}, {1.0, 1.0}, {1.0}), rl::InputError);
    CHECK_THROWS_AS(Tridiagonal<double>({}, {}, {}), rl::InputError);
    Tridiagonal<double> A({-1.0}, {2.0, 2.0}, {-1.0});
    CHECK_THROWS_AS(A.solve({1.0, 2.0, 3.0}), rl::InputError);
    CHECK_THROWS_AS(A.multiply({1.0}), rl::InputError);
}

TEST_CASE("Real solve equals double solve bit for bit", "[math][tridiag][aad]") {
    std::vector<double> a{-0.3, -0.4, -0.2}, b{2.0, 2.1, 2.2, 2.3}, c{-0.5, -0.1, -0.6}, d{1.0, -2.0, 3.0, 0.5};
    Tridiagonal<double> Ad(a, b, c);
    std::vector<Real> ar, br, cr, dr;
    for (double v : a) ar.push_back(Real::input(v));
    for (double v : b) br.push_back(Real::input(v));
    for (double v : c) cr.push_back(Real::input(v));
    for (double v : d) dr.push_back(Real::input(v));
    Tridiagonal<Real> Ar(ar, br, cr);
    auto xd = Ad.solve(d);
    auto xr = Ar.solve(dr);
    for (std::size_t i = 0; i < 4; ++i) {
        CHECK(xd[i] == xr[i].value());
        CHECK(xr[i].isActive());
    }
    rl::aad::Tape::current().reset();
}
