// SABR node calibration: the problem (residuals, Jacobian, bounds), the weighting, the starting point,
// recovery of known parameters, failure detection, and the fit of the shipped USD cube.
#include "rl/vol/sabr_fit.hpp"

#include "rl/core/error.hpp"
#include "rl/math/calibration_problem.hpp"
#include "rl/vol/cube.hpp"
#include "rl/vol/sabr.hpp"
#include "test_support.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <algorithm>
#include <cmath>

using namespace rl::vol;
using rl::market::SabrParams;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

static_assert(rl::math::CalibrationProblem<SabrNodeProblem>);

namespace {
std::vector<SmileQuote> syntheticQuotes(const test::NodeSpec& n, double roundBp = 0.0) {
    const auto p = n.params();
    std::vector<SmileQuote> q;
    for (double off : test::kOffsets) {
        double v = sabrNormalVol(p, n.forward, n.forward + off, n.expiry);
        if (off == 0.0) v = n.atm;
        if (roundBp > 0.0) v = std::round(v * 1e4 / roundBp) * roundBp / 1e4;
        q.push_back({n.forward + off, v});
    }
    return q;
}

SabrNodeProblem problemFor(const test::NodeSpec& n, SabrFitOptions o = {}, double roundBp = 0.0) {
    o.beta = n.beta;
    o.shift = n.shift;
    const auto q = syntheticQuotes(n, roundBp);
    const double atm = roundBp > 0.0 ? std::round(n.atm * 1e4 / roundBp) * roundBp / 1e4 : n.atm;
    return SabrNodeProblem("test", q, n.forward, n.expiry, atm, o);
}
} // namespace

// --- the problem ----------------------------------------------------------------------------------

TEST_CASE("SabrNodeProblem: layout, ids, sorting, and validation", "[vol][fit][problem]") {
    auto q = syntheticQuotes(test::kBenign);
    std::reverse(q.begin(), q.end());
    const SabrNodeProblem p("1Yx10Y", q, test::kBenign.forward, 1.0, test::kBenign.atm);
    CHECK(p.parameterCount() == 2);
    CHECK(p.residualCount() == 7);
    CHECK(p.quotes().front().strike < p.quotes().back().strike);
    REQUIRE(p.instrumentIds().size() == 7);
    CHECK(p.instrumentIds()[0] == "-100bp");
    CHECK(p.instrumentIds()[3] == "ATM");
    CHECK(p.instrumentIds()[4] == "+25bp");

    SabrFitOptions three;
    three.alphaFromAtm = false;
    CHECK(SabrNodeProblem("x", q, test::kBenign.forward, 1.0, test::kBenign.atm, three).parameterCount() == 3);

    const auto [lo, hi] = p.bounds();
    CHECK(lo.size() == 2);
    CHECK(lo[0] == p.options().rhoMin);
    CHECK(hi[1] == p.options().nuMax);

    CHECK_THROWS_AS(SabrNodeProblem("x", {q[0]}, 0.036, 1.0, 0.0088), rl::InputError);                                     // too few quotes
    CHECK_THROWS_AS(SabrNodeProblem("x", {q[0], q[0], q[1]}, 0.036, 1.0, 0.0088), rl::InputError);                         // duplicate strike
    CHECK_THROWS_AS(SabrNodeProblem("x", {{-0.001, 0.009}, {0.01, 0.009}}, 0.036, 1.0, 0.0088), rl::InputError);           // K + s ≤ 0, β > 0
    CHECK_THROWS_AS(SabrNodeProblem("x", q, 0.036, 0.0, 0.0088), rl::InputError);                                          // expiry
    SabrFitOptions bad;
    bad.rhoMin = 0.5;
    bad.rhoMax = 0.4;
    CHECK_THROWS_AS(SabrNodeProblem("x", q, 0.036, 1.0, 0.0088, bad), rl::InputError);
    // β = 0 accepts negative strikes.
    SabrFitOptions normal;
    normal.beta = 0.0;
    CHECK_NOTHROW(SabrNodeProblem("x", {{-0.001, 0.009}, {0.01, 0.009}}, 0.004, 1.0, 0.0088, normal));
}

TEST_CASE("SabrNodeProblem: residuals vanish at the generating parameters; params and toVector round-trip", "[vol][fit][problem]") {
    for (const auto& node : {test::kBenign, test::kStressed}) {
        const auto prob = problemFor(node);
        const auto z = prob.toVector(node.params());
        const auto r = prob.residuals(z);
        CHECK(r.lpNorm<Eigen::Infinity>() < 1e-10);
        const auto back = prob.params(z);
        CHECK_THAT(back.alpha, WithinRel(node.params().alpha, 1e-15));
        CHECK(back.rho == node.rho);
        CHECK_THROWS_AS(prob.params(rl::math::Vector::Zero(3)), rl::InputError);
    }
}

TEST_CASE("SabrNodeProblem: analytic Jacobian matches finite differences, with and without the ATM constraint", "[vol][fit][problem]") {
    for (bool fromAtm : {true, false}) {
        SabrFitOptions o;
        o.alphaFromAtm = fromAtm;
        const auto prob = problemFor(test::kBenign, o);
        rl::math::Vector z = prob.toVector(test::kBenign.params());
        z[z.size() - 1] *= 1.1; // off the solution, so the Jacobian is generic
        const auto J = prob.jacobian(z);
        const auto Jfd = rl::math::finiteDifferenceJacobian([&](const rl::math::Vector& x) { return prob.residuals(x); }, z, 1e-6);
        REQUIRE(J.rows() == Jfd.rows());
        REQUIRE(J.cols() == Jfd.cols());
        for (Eigen::Index i = 0; i < J.rows(); ++i)
            for (Eigen::Index j = 0; j < J.cols(); ++j) CHECK(std::fabs(J(i, j) - Jfd(i, j)) <= 1e-5 * (1.0 + std::fabs(Jfd(i, j))));
    }
}

TEST_CASE("smileQuotes: the ATM point is inserted and quotes are sorted by strike", "[vol][fit][problem]") {
    rl::market::SwaptionVolNodeQuote n;
    n.expiry = rl::core::Tenor::years(1);
    n.tenor = rl::core::Tenor::years(10);
    n.forward = 0.036;
    n.atmNormalVol = 0.0088;
    n.smile = {{-0.005, 0.0091}, {0.005, 0.0087}};
    const auto q = smileQuotes(n);
    REQUIRE(q.size() == 3);
    CHECK(q[1].strike == 0.036);
    CHECK(q[1].normalVol == 0.0088);
    CHECK_THAT(q[0].strike, WithinAbs(0.031, 1e-16));
    const auto p = makeNodeProblem(n, SabrFitOptions{});
    CHECK(p.id() == "1Yx10Y");
    CHECK(p.expiry() == 1.0);
}

// --- weighting, initial guess, fit -------------------------------------------------------------------

TEST_CASE("Weights: finite and positive; Equal is ones; Vega peaks at the money", "[vol][fit][weights]") {
    for (auto w : {SmileWeighting::Equal, SmileWeighting::Vega, SmileWeighting::InverseVol}) {
        SabrFitOptions o;
        o.weighting = w;
        const auto prob = problemFor(test::kBenign, o);
        const auto s = prob.scaling(prob.toVector(test::kBenign.params()));
        REQUIRE(static_cast<std::size_t>(s.size()) == prob.residualCount());
        for (Eigen::Index i = 0; i < s.size(); ++i) {
            CHECK(std::isfinite(s[i]));
            CHECK(s[i] > 0.0);
        }
        if (w == SmileWeighting::Equal) CHECK((s.array() == 1.0).all());
        if (w == SmileWeighting::Vega) {
            CHECK(s[3] >= s[2]);
            CHECK(s[3] >= s[4]);
            CHECK(s[2] >= s[0]);
            CHECK(s[4] >= s[6]);
        }
    }
}

TEST_CASE("Initial guess: inside the bounds and in the right neighbourhood", "[vol][fit][guess]") {
    for (const auto& node : {test::kBenign, test::kStressed}) {
        const auto prob = problemFor(node);
        const auto z = sabrInitialGuess(prob);
        const auto [lo, hi] = prob.bounds();
        REQUIRE(z.size() == 2);
        CHECK(z.allFinite());
        CHECK((z.array() > lo.array()).all());
        CHECK((z.array() < hi.array()).all());
        CHECK_THAT(z[0], WithinAbs(node.rho, 0.2));
        CHECK_THAT(z[1], WithinRel(node.nu, 0.5));
    }
}

TEST_CASE("Fit recovers the generating parameters of a synthetic smile", "[vol][fit]") {
    for (const auto& node : {test::kBenign, test::kStressed}) {
        for (bool fromAtm : {true, false}) {
            SabrFitOptions o;
            o.alphaFromAtm = fromAtm;
            const auto res = fitSabrSmile(problemFor(node, o));
            INFO("rho " << node.rho << " nu " << node.nu << " alphaFromAtm " << fromAtm);
            CHECK(res.ok);
            CHECK(res.flags.empty());
            CHECK(res.maxAbsBp < 1e-8);
            CHECK_THAT(res.params.rho, WithinAbs(node.rho, 1e-6));
            CHECK_THAT(res.params.nu, WithinAbs(node.nu, 1e-6));
            CHECK_THAT(res.params.alpha, WithinRel(node.params().alpha, 1e-6));
            CHECK(res.residualsBp.size() == 7);
            CHECK(res.status != rl::math::LmStatus::MaxIterations);
        }
    }
}

TEST_CASE("Fit to a smile quoted to 0.01 bp: under 0.1 bp everywhere and ATM repriced exactly", "[vol][fit]") {
    const auto prob = problemFor(test::kBenign, SabrFitOptions{}, 0.01);
    const auto res = fitSabrSmile(prob);
    CHECK(res.ok);
    CHECK(res.maxAbsBp < 0.1);
    CHECK(res.rmsBp <= res.maxAbsBp);
    CHECK_THAT(sabrNormalVol(res.params, prob.forward(), prob.forward(), prob.expiry()) * 1e4, WithinAbs(prob.atmNormalVol() * 1e4, 1e-10));
}

TEST_CASE("Fit failure detection: a frown is not a SABR smile", "[vol][fit][failure]") {
    // Wings BELOW the ATM vol: concave in strike, which no ν ≥ 0 produces.
    std::vector<SmileQuote> frown{{0.026, 0.0070}, {0.031, 0.0082}, {0.036, 0.0088}, {0.041, 0.0082}, {0.046, 0.0070}};
    SabrFitOptions o;
    const SabrNodeProblem prob("frown", frown, 0.036, 1.0, 0.0088, o);
    const auto res = fitSabrSmile(prob);
    CHECK_FALSE(res.ok);
    CHECK_FALSE(res.flags.empty());
    CHECK(res.maxAbsBp > o.maxResidualBp);
    o.throwOnFailure = true;
    CHECK_THROWS_AS(fitSabrSmile(SabrNodeProblem("frown", frown, 0.036, 1.0, 0.0088, o)), rl::CalibrationError);
}

TEST_CASE("Fit the shipped USD cube: every node under 0.1 bp, and the cube reprices each ATM exactly", "[vol][fit][data]") {
    const auto quotes = test::usdVolQuotes();
    const CubeFit fit = fitCube(quotes, SabrFitOptions{});
    REQUIRE(fit.fits.size() == 81);
    double worst = 0.0;
    for (std::size_t i = 0; i < fit.fits.size(); ++i) {
        INFO(fit.nodeIds[i]);
        CHECK(fit.fits[i].ok);
        worst = std::max(worst, fit.fits[i].maxAbsBp);
    }
    CHECK(worst < 0.1);
    for (std::size_t i = 0; i < fit.cube.expiryCount(); ++i) {
        for (std::size_t j = 0; j < fit.cube.tenorCount(); ++j) {
            const auto& n = fit.cube.node(i, j);
            const double e = fit.cube.expiries()[i], t = fit.cube.tenors()[j];
            CHECK_THAT(cubeNormalVol(fit.cube, e, t, n.forward, n.forward), WithinRel(n.atmNormalVol, 1e-12));
        }
    }
}
