#include "rl/aad/eigen_traits.hpp"

#include <catch2/catch_test_macros.hpp>

using rl::aad::MatrixR;
using rl::aad::Real;
using rl::aad::Tape;
using rl::aad::VectorR;

TEST_CASE("Eigen products over Real match double", "[aad][eigen]") {
    Tape& t = Tape::current();
    t.reset();
    Eigen::MatrixXd A(3, 3);
    A << 1, 2, 3, 4, 5, 6, 7, 8, 10;
    Eigen::VectorXd v(3);
    v << 0.5, -1.0, 2.0;

    MatrixR Ar(3, 3);
    VectorR vr(3);
    for (int i = 0; i < 3; ++i) {
        vr(i) = Real::input(v(i));
        for (int j = 0; j < 3; ++j) Ar(i, j) = Real(A(i, j));
    }
    VectorR y = Ar * vr;
    Eigen::VectorXd yd = A * v;
    for (int i = 0; i < 3; ++i) {
        CHECK(y(i).value() == yd(i));
        CHECK(y(i).isActive());
    }
    Real dot = vr.dot(vr);
    CHECK(dot.value() == v.dot(v));
    MatrixR AtA = Ar.transpose() * Ar;
    CHECK(AtA(2, 2).value() == (A.transpose() * A)(2, 2));
    CHECK_FALSE(AtA(2, 2).isActive());
    t.reset();
}
