// SwaptionCube<T>: validation, bilinear read with total-variance ATM, flat extrapolation, JSON, hashing.
#include "rl/market/swaption_cube.hpp"

#include "rl/aad/real.hpp"
#include "rl/core/error.hpp"
#include "rl/market/hash.hpp"

#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <cmath>
#include <string>
#include <vector>

using namespace rl::market;
using rl::aad::Real;
using rl::core::Date;
using rl::core::Tenor;
using Catch::Matchers::ContainsSubstring;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace {
const Date kAsOf(2026, 9, 4);

const double kAtm[3][2] = {{0.0100, 0.0092}, {0.0096, 0.0088}, {0.0084, 0.0078}};
const double kRho[3][2] = {{-0.05, -0.10}, {-0.10, -0.15}, {-0.20, -0.25}};
const double kNu[3][2] = {{0.60, 0.55}, {0.45, 0.40}, {0.25, 0.20}};

template <class T>
CubeNode<T> makeNode(double fwd, double atm, double alpha, double rho, double nu, double beta = 0.5) {
    CubeNode<T> n;
    n.forward = fwd;
    n.atmNormalVol = T(atm);
    n.sabr.alpha = T(alpha);
    n.sabr.rho = T(rho);
    n.sabr.nu = T(nu);
    n.sabr.beta = beta;
    n.sabr.shift = 0.0;
    return n;
}

std::vector<Tenor> expiries() { return {Tenor::months(6), Tenor::years(1), Tenor::years(5)}; }
std::vector<Tenor> tenors() { return {Tenor::years(2), Tenor::years(10)}; }

template <class T>
std::vector<CubeNode<T>> nodes() {
    std::vector<CubeNode<T>> out;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 2; ++j) out.push_back(makeNode<T>(0.035, kAtm[i][j], 0.05, kRho[i][j], kNu[i][j]));
    return out;
}

template <class T>
SwaptionCube<T> smallCube() { return SwaptionCube<T>(kAsOf, "USD-SOFR", expiries(), tenors(), nodes<T>(), 0.5, 0.0); }
} // namespace

TEST_CASE("Cube axis years and node labels", "[market][cube]") {
    CHECK(cubeAxisYears(Tenor::months(1)) == 1.0 / 12.0);
    CHECK(cubeAxisYears(Tenor::months(6)) == 0.5);
    CHECK(cubeAxisYears(Tenor::years(10)) == 10.0);
    CHECK(cubeAxisYears(Tenor::weeks(1)) == 7.0 / 365.0);
    CHECK(cubeAxisYears(Tenor::days(1)) == 1.0 / 365.0);
    CHECK_THROWS_AS(cubeAxisYears(Tenor::months(0)), rl::InputError);
    CHECK(cubeNodeLabel(Tenor::years(1), Tenor::years(10)) == "1Yx10Y");
}

TEST_CASE("Cube validation rejects malformed grids and parameters", "[market][cube]") {
    const auto good = smallCube<double>();
    CHECK(good.size() == 6);
    CHECK(good.expiryCount() == 3);
    CHECK(good.tenorCount() == 2);
    CHECK(good.expiries()[0] == 0.5);

    auto build = [](std::vector<Tenor> ex, std::vector<Tenor> tn, std::vector<CubeNode<double>> ns, double beta = 0.5) {
        return SwaptionCube<double>(kAsOf, "USD-SOFR", std::move(ex), std::move(tn), std::move(ns), beta, 0.0);
    };
    CHECK_THROWS_AS(build({Tenor::years(1)}, tenors(), {nodes<double>()[0], nodes<double>()[1]}), rl::InputError);
    CHECK_THROWS_AS(build({Tenor::years(1), Tenor::months(6), Tenor::years(5)}, tenors(), nodes<double>()), rl::InputError);
    CHECK_THROWS_AS(build({Tenor::months(12), Tenor::years(1), Tenor::years(5)}, tenors(), nodes<double>()), rl::InputError); // same point
    auto five = nodes<double>();
    five.pop_back();
    CHECK_THROWS_AS(build(expiries(), tenors(), five), rl::InputError);

    auto mutate = [&](auto f) {
        auto ns = nodes<double>();
        f(ns[3]);
        return build(expiries(), tenors(), ns);
    };
    CHECK_THROWS_AS(mutate([](CubeNode<double>& n) { n.sabr.rho = 1.0; }), rl::InputError);
    CHECK_THROWS_AS(mutate([](CubeNode<double>& n) { n.sabr.alpha = 0.0; }), rl::InputError);
    CHECK_THROWS_AS(mutate([](CubeNode<double>& n) { n.sabr.nu = -0.1; }), rl::InputError);
    CHECK_THROWS_AS(mutate([](CubeNode<double>& n) { n.sabr.beta = 0.3; }), rl::InputError);
    CHECK_THROWS_AS(mutate([](CubeNode<double>& n) { n.atmNormalVol = 0.0; }), rl::InputError);
    CHECK_THROWS_AS(mutate([](CubeNode<double>& n) { n.forward = -0.01; }), rl::InputError); // β > 0 needs F + s > 0
    try {
        mutate([](CubeNode<double>& n) { n.sabr.rho = -1.5; });
        FAIL("expected InputError");
    } catch (const rl::InputError& e) {
        CHECK_THAT(std::string(e.what()), ContainsSubstring("1Yx10Y")); // node 3 is (1Y, 10Y): errors name the node
    }
    // β = 0 (normal SABR) allows a negative forward.
    auto neg = nodes<double>();
    for (auto& n : neg) n.sabr.beta = 0.0;
    neg[0].forward = -0.005;
    CHECK_NOTHROW(build(expiries(), tenors(), neg, 0.0));
}

TEMPLATE_TEST_CASE("Cube interpolation: nodes reproduced, total-variance ATM, bilinear smile, flat beyond the axes", "[market][cube]", double, Real) {
    using T = TestType;
    const auto c = smallCube<T>();
    for (std::size_t i = 0; i < 3; ++i) {
        for (std::size_t j = 0; j < 2; ++j) {
            const auto in = c.interpolate(c.expiries()[i], c.tenors()[j]);
            CHECK_THAT(rl::math::value(in.atmNormalVol), WithinRel(kAtm[i][j], 1e-14));
            CHECK_THAT(rl::math::value(in.rho), WithinAbs(kRho[i][j], 1e-15));
            CHECK_THAT(rl::math::value(in.nu), WithinAbs(kNu[i][j], 1e-15));
            CHECK(in.beta == 0.5);
            CHECK(in.shift == 0.0);
        }
    }
    // Halfway in expiry between 1Y and 5Y at the 2Y tenor: total variance is linear in expiry.
    {
        const auto in = c.interpolate(3.0, 2.0);
        const double var = 0.5 * kAtm[1][0] * kAtm[1][0] * 1.0 + 0.5 * kAtm[2][0] * kAtm[2][0] * 5.0;
        CHECK_THAT(rl::math::value(in.atmNormalVol), WithinRel(std::sqrt(var / 3.0), 1e-14));
        CHECK_THAT(rl::math::value(in.rho), WithinAbs(-0.15, 1e-15));
        // Between the neighbouring nodes' vols, as a total-variance interpolation must be.
        CHECK(rl::math::value(in.atmNormalVol) < kAtm[1][0]);
        CHECK(rl::math::value(in.atmNormalVol) > kAtm[2][0]);
    }
    // Centre of the (1Y..5Y) × (2Y..10Y) cell: every weight is ¼.
    {
        const auto in = c.interpolate(3.0, 6.0);
        CHECK_THAT(rl::math::value(in.nu), WithinAbs(0.25 * (0.45 + 0.40 + 0.25 + 0.20), 1e-15));
        CHECK_THAT(rl::math::value(in.rho), WithinAbs(0.25 * (-0.10 - 0.15 - 0.20 - 0.25), 1e-15));
    }
    // Flat extrapolation in both directions and both axes.
    {
        const auto lo = c.interpolate(0.1, 30.0);
        CHECK_THAT(rl::math::value(lo.atmNormalVol), WithinRel(kAtm[0][1], 1e-14));
        CHECK_THAT(rl::math::value(lo.nu), WithinAbs(kNu[0][1], 1e-15));
        const auto hi = c.interpolate(20.0, 1.0);
        CHECK_THAT(rl::math::value(hi.atmNormalVol), WithinRel(kAtm[2][0], 1e-14));
        CHECK_THAT(rl::math::value(hi.rho), WithinAbs(kRho[2][0], 1e-15));
    }
    CHECK_THROWS_AS(c.interpolate(0.0, 5.0), rl::InputError);
    CHECK_THROWS_AS(c.interpolate(1.0, -1.0), rl::InputError);
    rl::aad::Tape::current().reset();
}

TEST_CASE("Cube interpolation under Real: active inputs give an active read", "[market][cube][aad]") {
    auto ns = nodes<Real>();
    for (auto& n : ns) n.atmNormalVol = Real::input(rl::math::value(n.atmNormalVol));
    const SwaptionCube<Real> c(kAsOf, "USD-SOFR", expiries(), tenors(), ns, 0.5, 0.0);
    const auto in = c.interpolate(2.0, 5.0);
    CHECK(in.atmNormalVol.isActive());
    CHECK(in.atmNormalVol.value() == smallCube<double>().interpolate(2.0, 5.0).atmNormalVol);
    rl::aad::Tape::current().reset();
}

TEST_CASE("Cube: node access, indexOf, withNode", "[market][cube]") {
    const auto c = smallCube<double>();
    CHECK(c.node(1, 1).atmNormalVol == kAtm[1][1]);
    CHECK_THROWS_AS(c.node(3, 0), rl::InputError);
    CHECK(c.indexOf(Tenor::years(1), Tenor::years(10)) == std::pair<std::size_t, std::size_t>{1, 1});
    CHECK(c.indexOf(Tenor::months(12), Tenor::years(10)) == std::pair<std::size_t, std::size_t>{1, 1}); // 12M == 1Y
    CHECK_THROWS_AS(c.indexOf(Tenor::years(2), Tenor::years(10)), rl::InputError);
    auto n = c.node(1, 1);
    n.sabr.nu = 0.9;
    const auto c2 = c.withNode(1, 1, n);
    CHECK(c2.node(1, 1).sabr.nu == 0.9);
    CHECK(c.node(1, 1).sabr.nu == kNu[1][1]); // original untouched
    n.sabr.rho = 2.0;
    CHECK_THROWS_AS(c.withNode(1, 1, n), rl::InputError); // re-validated
}

TEST_CASE("Cube JSON: round trip, any node order, loud on holes and duplicates", "[market][cube][json]") {
    const auto c = smallCube<double>();
    const auto back = swaptionCubeFromJson(toJson(c));
    CHECK(back.asOf() == kAsOf);
    CHECK(back.index() == "USD-SOFR");
    CHECK(back.beta() == 0.5);
    REQUIRE(back.size() == c.size());
    for (std::size_t k = 0; k < c.size(); ++k) {
        CHECK_THAT(back.nodes()[k].atmNormalVol, WithinRel(c.nodes()[k].atmNormalVol, 1e-15));
        CHECK_THAT(back.nodes()[k].forward, WithinRel(c.nodes()[k].forward, 1e-15));
        CHECK(back.nodes()[k].sabr.rho == c.nodes()[k].sabr.rho);
        CHECK(back.nodes()[k].sabr.nu == c.nodes()[k].sabr.nu);
    }
    const std::string hdr = R"({"asOf": "2026-09-04", "index": "USD-SOFR", "beta": 0.5, "shift": 0.0, "expiries": ["1Y", "5Y"], "tenors": ["2Y", "10Y"], "nodes": [)";
    auto node = [](const char* e, const char* t) {
        return std::string(R"({"expiry": ")") + e + R"(", "tenor": ")" + t + R"(", "forward_pct": 3.5, "atm_bp": 90, "alpha": 0.05, "rho": -0.1, "nu": 0.4})";
    };
    const std::string shuffled = hdr + node("5Y", "10Y") + "," + node("1Y", "2Y") + "," + node("5Y", "2Y") + "," + node("1Y", "10Y") + "]}";
    const auto s = swaptionCubeFromJson(shuffled);
    CHECK(s.node(0, 1).atmNormalVol == 0.009);
    CHECK_THROWS_AS(swaptionCubeFromJson(hdr + node("5Y", "10Y") + "," + node("1Y", "2Y") + "," + node("5Y", "2Y") + "]}"), rl::InputError);
    CHECK_THROWS_AS(swaptionCubeFromJson(hdr + node("5Y", "10Y") + "," + node("5Y", "10Y") + "," + node("1Y", "2Y") + "," + node("5Y", "2Y") + "," +
                                         node("1Y", "10Y") + "]}"),
                    rl::InputError);
    CHECK_THROWS_AS(swaptionCubeFromJson(hdr + node("3Y", "10Y") + "]}"), rl::InputError);
    CHECK_THROWS_AS(swaptionCubeFromJson("{not json"), rl::InputError);
}

TEST_CASE("Cube hash: content-determined, one-ulp sensitive, the same under Real", "[market][cube][hash]") {
    auto hashOf = [](const auto& cube) {
        Hasher h;
        cube.mixInto(h);
        return h.value();
    };
    const auto base = hashOf(smallCube<double>());
    CHECK(hashOf(smallCube<double>()) == base);
    auto n = smallCube<double>().node(2, 1);
    n.sabr.nu = std::nextafter(n.sabr.nu, 1.0);
    CHECK(hashOf(smallCube<double>().withNode(2, 1, n)) != base);
    auto m = smallCube<double>().node(0, 0);
    m.atmNormalVol = std::nextafter(m.atmNormalVol, 1.0);
    CHECK(hashOf(smallCube<double>().withNode(0, 0, m)) != base);
    CHECK(hashOf(smallCube<Real>()) == base);
    rl::aad::Tape::current().reset();
}
