#include "rl/math/brownian_bridge.hpp"

#include "rl/core/error.hpp"
#include "rl/math/rng.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <Eigen/Dense>

#include <cmath>
#include <vector>

using rl::math::BrownianBridge;
using Catch::Matchers::WithinAbs;

TEST_CASE("Single step is the terminal value", "[math][bridge]") {
    BrownianBridge bb({2.25});
    CHECK(bb.size() == 1);
    CHECK_THAT(bb.stdDev()[0], WithinAbs(1.5, 1e-15));
    CHECK(bb.bridgeIndex()[0] == 0);
    std::vector<double> g{0.8}, w(1);
    bb.buildPath(g, w);
    CHECK_THAT(w[0], WithinAbs(1.2, 1e-15));
}

TEST_CASE("Uniform 8-step bridge: order and weights", "[math][bridge]") {
    std::vector<double> t;
    for (int i = 1; i <= 8; ++i) t.push_back(0.125 * i);
    BrownianBridge bb(t);
    CHECK(bb.bridgeIndex() == std::vector<std::size_t>{7, 3, 1, 5, 0, 2, 4, 6});
    // Step 1 fills t = 0.5 between 0 and t = 1: weights ½, ½, stdDev √(0.5·0.5/1) = 0.5.
    CHECK_THAT(bb.leftWeight()[1], WithinAbs(0.5, 1e-15));
    CHECK_THAT(bb.rightWeight()[1], WithinAbs(0.5, 1e-15));
    CHECK_THAT(bb.stdDev()[1], WithinAbs(0.5, 1e-15));
    CHECK(bb.rightIndex()[1] == 7);
    // Step 2 fills t = 0.25 between 0 and 0.5: stdDev √(0.25·0.25/0.5) = √0.125.
    CHECK_THAT(bb.stdDev()[2], WithinAbs(std::sqrt(0.125), 1e-15));
    CHECK(bb.rightIndex()[2] == 3);
    CHECK_THAT(bb.stdDev()[0], WithinAbs(1.0, 1e-15));
}

TEST_CASE("Covariance of the linear map is min(t_i, t_j) exactly", "[math][bridge]") {
    std::vector<double> t{0.1, 0.25, 0.5, 0.6, 1.0, 1.5, 2.0};
    const std::size_t n = t.size();
    BrownianBridge bb(t);
    Eigen::MatrixXd L(n, n);
    std::vector<double> g(n, 0.0), w(n);
    for (std::size_t k = 0; k < n; ++k) {
        std::fill(g.begin(), g.end(), 0.0);
        g[k] = 1.0;
        bb.buildPath(g, w);
        for (std::size_t i = 0; i < n; ++i) L(static_cast<int>(i), static_cast<int>(k)) = w[i];
    }
    Eigen::MatrixXd C = L * L.transpose();
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) CHECK_THAT(C(static_cast<int>(i), static_cast<int>(j)), WithinAbs(std::min(t[i], t[j]), 1e-12));
}

TEST_CASE("Increments are differences of the path", "[math][bridge]") {
    std::vector<double> t{0.5, 1.0, 1.5, 2.0, 2.5};
    BrownianBridge bb(t);
    rl::math::Rng rng(3);
    std::vector<double> g(5), w(5), dw(5);
    for (auto& x : g) x = rng.gaussian();
    bb.buildPath(g, w);
    bb.buildIncrements(g, dw);
    CHECK_THAT(dw[0], WithinAbs(w[0], 1e-15));
    for (std::size_t i = 1; i < 5; ++i) CHECK_THAT(dw[i], WithinAbs(w[i] - w[i - 1], 1e-15));
}

TEST_CASE("Statistical sanity: variance of increments equals dt", "[math][bridge]") {
    std::vector<double> t{0.25, 0.5, 1.0, 2.0};
    BrownianBridge bb(t);
    rl::math::Rng rng(11);
    const int N = 100000;
    std::vector<double> g(4), dw(4), s2(4, 0.0);
    for (int k = 0; k < N; ++k) {
        for (auto& x : g) x = rng.gaussian();
        bb.buildIncrements(g, dw);
        for (std::size_t i = 0; i < 4; ++i) s2[i] += dw[i] * dw[i];
    }
    CHECK_THAT(s2[0] / N, WithinAbs(0.25, 0.01));
    CHECK_THAT(s2[1] / N, WithinAbs(0.25, 0.01));
    CHECK_THAT(s2[2] / N, WithinAbs(0.5, 0.02));
    CHECK_THAT(s2[3] / N, WithinAbs(1.0, 0.03));
}

TEST_CASE("Validation", "[math][bridge]") {
    CHECK_THROWS_AS(BrownianBridge({}), rl::InputError);
    CHECK_THROWS_AS(BrownianBridge({0.0, 1.0}), rl::InputError);
    CHECK_THROWS_AS(BrownianBridge({1.0, 0.5}), rl::InputError);
    BrownianBridge bb({1.0, 2.0});
    std::vector<double> g(1), w(2);
    CHECK_THROWS_AS(bb.buildPath(g, w), rl::InputError);
}
