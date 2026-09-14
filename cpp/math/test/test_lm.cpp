#include "rl/math/lm.hpp"

#include "rl/math/rng.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using namespace rl::math;
using Catch::Matchers::WithinAbs;

TEST_CASE("Rosenbrock from the classic start", "[math][lm]") {
    ResidualFn r = [](const Vector& x) {
        Vector out(2);
        out << 10.0 * (x(1) - x(0) * x(0)), 1.0 - x(0);
        return out;
    };
    JacobianFn j = [](const Vector& x) {
        Matrix J(2, 2);
        J << -20.0 * x(0), 10.0, -1.0, 0.0;
        return J;
    };
    Vector x0(2);
    x0 << -1.2, 1.0;
    auto res = levenbergMarquardt(r, j, x0);
    CHECK(res.status != LmStatus::MaxIterations);
    CHECK(res.status != LmStatus::Failed);
    CHECK_THAT(res.x(0), WithinAbs(1.0, 1e-8));
    CHECK_THAT(res.x(1), WithinAbs(1.0, 1e-8));
    CHECK(res.cost < 1e-16);
    CHECK(res.iterations < 100);

    // Finite-difference Jacobian path gives the same answer.
    auto fd = levenbergMarquardt(r, {}, x0);
    CHECK_THAT(fd.x(0), WithinAbs(1.0, 1e-7));
    CHECK_THAT(fd.x(1), WithinAbs(1.0, 1e-7));
}

TEST_CASE("Linear least squares matches QR", "[math][lm]") {
    Rng rng(7);
    const int m = 12, n = 3;
    Matrix A(m, n);
    Vector b(m);
    for (int i = 0; i < m; ++i) {
        for (int k = 0; k < n; ++k) A(i, k) = rng.uniform01() - 0.5;
        b(i) = rng.uniform01();
    }
    ResidualFn r = [&](const Vector& x) { return Vector(A * x - b); };
    JacobianFn j = [&](const Vector&) { return A; };
    Vector ref = A.colPivHouseholderQr().solve(b);
    auto res = levenbergMarquardt(r, j, Vector::Zero(n));
    for (int k = 0; k < n; ++k) CHECK_THAT(res.x(k), WithinAbs(ref(k), 1e-10));
}

TEST_CASE("Gradient already zero returns immediately", "[math][lm]") {
    ResidualFn r = [](const Vector& x) { return Vector(x); };
    JacobianFn j = [](const Vector& x) { return Matrix(Matrix::Identity(x.size(), x.size())); };
    auto res = levenbergMarquardt(r, j, Vector::Zero(2));
    CHECK(res.iterations == 0);
    CHECK(res.status == LmStatus::GradientSmall);
}

TEST_CASE("Box bounds are respected", "[math][lm]") {
    // Minimum of (x - 2)² is at 2; bound x ≤ 1.
    ResidualFn r = [](const Vector& x) {
        Vector out(1);
        out << x(0) - 2.0;
        return out;
    };
    LmOptions opt;
    opt.lowerBounds = Vector::Constant(1, -10.0);
    opt.upperBounds = Vector::Constant(1, 1.0);
    auto res = levenbergMarquardt(r, {}, Vector::Zero(1), opt);
    CHECK_THAT(res.x(0), WithinAbs(1.0, 1e-9));
}

TEST_CASE("Finite-difference Jacobian on a quadratic", "[math][lm]") {
    Matrix A(3, 2);
    A << 1, 2, 3, 4, 5, 6;
    ResidualFn r = [&](const Vector& x) { return Vector(A * x + Vector::Constant(3, 0.5 * x.squaredNorm())); };
    Vector x(2);
    x << 0.3, -0.7;
    Matrix J = finiteDifferenceJacobian(r, x);
    // d/dx_k of 0.5 Σ x_i² = x_k, added to every row.
    for (int i = 0; i < 3; ++i)
        for (int k = 0; k < 2; ++k) CHECK_THAT(J(i, k), WithinAbs(A(i, k) + x(k), 1e-8));
}

TEST_CASE("Status names", "[math][lm]") {
    CHECK(std::string(toString(LmStatus::Converged)) == "Converged");
    CHECK(std::string(toString(LmStatus::MaxIterations)) == "MaxIterations");
}
