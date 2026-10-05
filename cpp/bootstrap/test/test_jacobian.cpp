// The adjoint Jacobian against finite differences, and par risk by the implicit-function theorem against
// bump-and-rebuild.
#include "rl/bootstrap/global_bootstrap.hpp"
#include "rl/bootstrap/jacobian.hpp"

#include "rl/instruments/builder.hpp"
#include "rl/linear/dv01.hpp"
#include "rl/linear/leg_pv.hpp"
#include "test_support.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>

using namespace rl::bootstrap;
using namespace rl::bootstrap::test;
using rl::core::Tenor;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace {
void checkClose(const rl::math::Matrix& a, const rl::math::Matrix& b, double rel, double absFloor, const BootstrapProblem& p) {
    REQUIRE(a.rows() == b.rows());
    REQUIRE(a.cols() == b.cols());
    for (Eigen::Index i = 0; i < a.rows(); ++i)
        for (Eigen::Index j = 0; j < a.cols(); ++j) {
            const double tol = absFloor + rel * std::max(std::fabs(a(i, j)), std::fabs(b(i, j)));
            INFO("row " << p.instrumentIds()[static_cast<std::size_t>(i)] << " col " << j);
            CHECK_THAT(a(i, j), WithinAbs(b(i, j), tol));
        }
}
} // namespace

TEST_CASE("dF/dz matches central differences at the guess and at the solution", "[bootstrap][jacobian]") {
    for (auto interp : {rl::market::CurveInterpolation::LinearLogDf, rl::market::CurveInterpolation::MonotoneConvex}) {
        BootstrapProblem p = usdProblem(interp);
        // At the flat initial guess every monotone-convex deviation g0 = g1 = 0 sits exactly on the
        // boundary between the four regimes, so the residual is only piecewise differentiable there
        // and a central difference straddles the kink. Compare at the guess for linear log-DF only.
        if (interp == rl::market::CurveInterpolation::LinearLogDf) {
            const rl::math::Vector z0 = initialGuess(p);
            checkClose(dFdz(p, z0), finiteDifferenceJacobian(p, z0), 1e-6, 1e-10, p);
        }
        BootstrapResult r = solve(p);
        REQUIRE(r.converged);
        // Monotone convex is C1 but not C2 in the nodes where a deviation pair sits on a regime
        // boundary, and the end intervals sit on one by construction (f_0 is defined so that
        // g0 = -g1/2). A central difference straddling such a point is only first-order accurate,
        // so the finite-difference reference is good to ~1e-3 relative there; the adjoint is exact.
        const bool linear = interp == rl::market::CurveInterpolation::LinearLogDf;
        checkClose(dFdz(p, r.z), finiteDifferenceJacobian(p, r.z), linear ? 1e-6 : 1e-3, linear ? 1e-10 : 1e-7, p);
    }
}

TEST_CASE("dF/dz structure: SOFR-only rows do not see FF nodes; basis rows see both curves", "[bootstrap][jacobian]") {
    BootstrapProblem p = usdProblem();
    BootstrapResult r = solve(p);
    REQUIRE(r.converged);
    const rl::math::Matrix J = dFdz(p, r.z);
    const auto [ffOff, ffN] = p.parameterRange("USD-FEDFUNDS");
    const auto [sofrOff, sofrN] = p.parameterRange("USD-SOFR");
    for (std::size_t i = 0; i < p.size(); ++i) {
        const auto& bi = p.instruments()[i];
        const auto row = J.row(static_cast<Eigen::Index>(i));
        const double ffPart = row.segment(static_cast<Eigen::Index>(ffOff), static_cast<Eigen::Index>(ffN)).cwiseAbs().maxCoeff();
        const double sofrPart = row.segment(static_cast<Eigen::Index>(sofrOff), static_cast<Eigen::Index>(sofrN)).cwiseAbs().maxCoeff();
        INFO(bi.id);
        if (bi.type == rl::market::QuoteType::BasisSwap) {
            CHECK(ffPart > 0.0);
            CHECK(sofrPart > 0.0);
        } else if (bi.curveName == "USD-SOFR") {
            CHECK(ffPart == 0.0);
        } else { // FF OIS: FF forwards and SOFR discounting
            CHECK(ffPart > 0.0);
            CHECK(sofrPart > 0.0);
        }
    }
    // Non-singular.
    Eigen::PartialPivLU<rl::math::Matrix> lu(J);
    CHECK(lu.rcond() > 1e-12);
}

TEST_CASE("Newton with the analytic Jacobian and with finite differences reach the same z", "[bootstrap][jacobian]") {
    BootstrapProblem p = usdProblem();
    BootstrapOptions a, f;
    a.analyticJacobian = true;
    f.analyticJacobian = false;
    BootstrapResult ra = solve(p, a), rf = solve(p, f);
    REQUIRE(ra.converged);
    REQUIRE(rf.converged);
    CHECK((ra.z - rf.z).cwiseAbs().maxCoeff() < 1e-12);
}

TEST_CASE("dF/dq matches central differences over withQuotes", "[bootstrap][jacobian]") {
    BootstrapProblem p = usdProblem();
    BootstrapResult r = solve(p);
    REQUIRE(r.converged);
    const rl::math::Matrix Jq = dFdq(p, r.z);
    const rl::math::Vector q = p.quoteValues();
    const double h = 1e-6;
    rl::math::Matrix fd(Jq.rows(), Jq.cols());
    for (Eigen::Index j = 0; j < q.size(); ++j) {
        rl::math::Vector qp = q, qm = q;
        qp[j] += h;
        qm[j] -= h;
        fd.col(j) = (p.withQuotes(qp).residuals(r.z) - p.withQuotes(qm).residuals(r.z)) / (2 * h);
    }
    checkClose(Jq, fd, 1e-7, 1e-9, p);
    // Diagonal: each quote is a contract term of exactly one instrument.
    for (Eigen::Index i = 0; i < Jq.rows(); ++i)
        for (Eigen::Index j = 0; j < Jq.cols(); ++j)
            if (i != j) CHECK(Jq(i, j) == 0.0);
}

TEST_CASE("Par risk via the implicit-function theorem matches bump-and-rebuild to 1e-6 of DV01", "[bootstrap][parrisk]") {
    BootstrapProblem p = usdProblem();
    BootstrapResult r = solve(p);
    REQUIRE(r.converged);
    // The trade: a 10y SOFR payer at 3.00% (off-market), notional 1e6.
    auto trade = rl::instruments::OisSwapBuilder::usdSofr(p.asOf()).id("T").tenor(Tenor::years(10)).fixedRate(0.03).notional(1e6).build();
    const auto& fixings = p.fixings();
    PvFunction pv = [&](const rl::market::CurveSet<double>& c) { return rl::linear::pv(trade, c, fixings); };
    const rl::math::Vector g = zeroRiskByBump(p, r.z, pv);
    const rl::math::Vector ift = parRisk(p, r.z, g);
    const rl::math::Vector bump = parRiskByRebuild(p, r.z, pv, 1e-6); // small bump: first order
    REQUIRE(ift.size() == bump.size());
    const double dv01 = std::fabs(rl::linear::fixedRateDv01(trade, p.curves(r.z))) * 1e4; // per unit of rate
    double total = 0.0;
    for (Eigen::Index i = 0; i < ift.size(); ++i) {
        INFO(p.instrumentIds()[static_cast<std::size_t>(i)]);
        CHECK_THAT(ift[i], WithinAbs(bump[i], 1e-6 * dv01));
        if (p.instruments()[static_cast<std::size_t>(i)].source.unit != rl::market::QuoteUnit::Price) total += ift[i];
    }
    // The 10y SOFR OIS carries most of the risk; the 25Y basis carries essentially none.
    const auto i10 = static_cast<Eigen::Index>(std::find(p.instrumentIds().begin(), p.instrumentIds().end(), "SOFR-OIS-10Y") - p.instrumentIds().begin());
    const auto i25 = static_cast<Eigen::Index>(std::find(p.instrumentIds().begin(), p.instrumentIds().end(), "FF-SOFR-25Y") - p.instrumentIds().begin());
    CHECK(std::fabs(ift[i10]) > 0.5 * dv01);
    CHECK(std::fabs(ift[i25]) < 1e-3 * dv01);
    // Sum of rate-quote par risks ≈ parallel DV01 in rate units (the futures are priced inversely, excluded above).
    CHECK_THAT(total, WithinRel(rl::linear::parallelDv01(trade, p.curves(r.z), fixings) * 1e4, 0.15));
}
