// Cube queries: node reproduction, ATM repriced off-grid at any forward, smile sharing, Real, and the snapshot.
#include "rl/vol/cube.hpp"

#include "rl/aad/real.hpp"
#include "rl/market/curve_set.hpp"
#include "rl/market/snapshot_builder.hpp"
#include "rl/market/yield_curve.hpp"
#include "rl/vol/sabr.hpp"

#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <vector>

using namespace rl::vol;
using rl::aad::Real;
using rl::core::Date;
using rl::core::Tenor;
using rl::market::CubeNode;
using rl::market::SabrParams;
using rl::market::SwaptionCube;
using rl::math::value;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace {
const Date kAsOf(2026, 9, 4);
const std::vector<Tenor> kEx{Tenor::months(6), Tenor::years(2), Tenor::years(10)};
const std::vector<Tenor> kTn{Tenor::years(2), Tenor::years(10), Tenor::years(30)};

struct Raw {
    double fwd, atm, rho, nu;
};
Raw raw(std::size_t i, std::size_t j) {
    const double e[] = {0.5, 2.0, 10.0}, t[] = {2.0, 10.0, 30.0};
    return {0.033 + 0.0015 * std::log1p(e[i] + t[j]) / 3.5, 0.0100 - 0.0010 * i - 0.0004 * j, -0.05 - 0.08 * i - 0.02 * j, 0.60 - 0.18 * i - 0.03 * j};
}

template <class T>
SwaptionCube<T> makeCube() {
    std::vector<CubeNode<T>> nodes;
    for (std::size_t i = 0; i < 3; ++i) {
        for (std::size_t j = 0; j < 3; ++j) {
            const Raw r = raw(i, j);
            CubeNode<T> n;
            n.forward = r.fwd;
            n.atmNormalVol = T(r.atm);
            n.sabr.rho = T(r.rho);
            n.sabr.nu = T(r.nu);
            n.sabr.beta = 0.5;
            n.sabr.shift = 0.0;
            n.sabr.alpha = alphaFromAtmNormalVol<T>(T(r.atm), T(r.fwd), rl::market::cubeAxisYears(kEx[i]), 0.5, T(r.rho), T(r.nu), 0.0);
            nodes.push_back(n);
        }
    }
    return SwaptionCube<T>(kAsOf, "USD-SOFR", kEx, kTn, nodes, 0.5, 0.0);
}
} // namespace

TEST_CASE("Cube query at a node reproduces the node's own SABR smile at its forward", "[vol][cube]") {
    const auto cube = makeCube<double>();
    for (std::size_t i = 0; i < 3; ++i) {
        for (std::size_t j = 0; j < 3; ++j) {
            const auto& n = cube.node(i, j);
            const double e = cube.expiries()[i], t = cube.tenors()[j];
            const auto p = cubeParams(cube, e, t, n.forward);
            CHECK_THAT(p.alpha, WithinRel(n.sabr.alpha, 1e-13));
            for (double off : {-0.01, -0.0025, 0.0, 0.005, 0.01}) {
                CHECK_THAT(cubeNormalVol(cube, e, t, n.forward + off, n.forward), WithinRel(sabrNormalVol(n.sabr, n.forward, n.forward + off, e), 1e-13));
            }
        }
    }
}

TEST_CASE("Cube reprices the interpolated ATM vol at any (expiry, tenor) and any forward", "[vol][cube]") {
    const auto cube = makeCube<double>();
    for (auto [e, t] : {std::pair{1.0, 5.0}, std::pair{3.7, 22.0}, std::pair{0.1, 40.0}, std::pair{15.0, 1.0}}) {
        const double atm = cubeAtmNormalVol(cube, e, t);
        for (double fwd : {0.020, 0.035, 0.050}) {
            INFO("e " << e << " t " << t << " F " << fwd);
            CHECK_THAT(cubeNormalVol(cube, e, t, fwd, fwd), WithinRel(atm, 1e-13));
        }
    }
    // The smile moves with the forward through α (the backbone): at a fixed strike the vol changes, ATM does not.
    CHECK(std::fabs(cubeNormalVol(cube, 5.0, 10.0, 0.035, 0.030) - cubeNormalVol(cube, 5.0, 10.0, 0.035, 0.035)) > 1e-6);
}

TEST_CASE("cubeSmile equals strike-by-strike queries", "[vol][cube]") {
    const auto cube = makeCube<double>();
    const std::vector<double> ks{0.02, 0.03, 0.04, 0.05};
    const auto smile = cubeSmile(cube, 4.0, 12.0, 0.036, ks);
    REQUIRE(smile.size() == ks.size());
    for (std::size_t i = 0; i < ks.size(); ++i) CHECK(smile[i] == cubeNormalVol(cube, 4.0, 12.0, ks[i], 0.036));
}

TEST_CASE("Cube query under Real: bit-identical and active in the forward and the ATM inputs", "[vol][cube][aad]") {
    const auto cd = makeCube<double>();
    auto cr = makeCube<Real>();
    auto n = cr.node(1, 1);
    n.atmNormalVol = Real::input(value(n.atmNormalVol));
    cr = cr.withNode(1, 1, n);
    const Real fwd = Real::input(0.036);
    const Real v = cubeNormalVol(cr, 3.0, 12.0, 0.03, fwd);
    CHECK(v.isActive());
    CHECK(v.value() == cubeNormalVol(cd, 3.0, 12.0, 0.03, 0.036));
    rl::aad::Tape::current().reset();
}

TEST_CASE("Cube carried in a snapshot is queryable through the snapshot", "[vol][cube][snapshot]") {
    std::vector<double> ts{1.0, 10.0, 30.0};
    std::vector<double> ys{0.035, 0.36, 1.12};
    rl::market::CurveSet<double> curves({rl::market::YieldCurve<double>(kAsOf, "USD-SOFR-OIS", ts, ys, rl::market::CurveInterpolation::LinearLogDf)},
                                        "USD-SOFR-OIS", {{"USD-SOFR", "USD-SOFR-OIS"}});
    const auto snap = rl::market::SnapshotBuilder<double>().asOf(kAsOf).curves(curves).swaptionCube(makeCube<double>()).build();
    REQUIRE(snap->swaptionCube() != nullptr);
    CHECK_THAT(cubeNormalVol(*snap->swaptionCube(), 2.0, 10.0, 0.035, 0.035), WithinRel(cubeAtmNormalVol(makeCube<double>(), 2.0, 10.0), 1e-13));
}
