// End-to-end checks of the global bootstrap on the shipped USD quote set.
#include "rl/bootstrap/global_bootstrap.hpp"
#include "rl/bootstrap/report.hpp"
#include "rl/bootstrap/residuals.hpp"

#include "rl/core/error.hpp"
#include "rl/linear/forward.hpp"
#include "rl/linear/leg_pv.hpp"
#include "rl/market/snapshot_builder.hpp"
#include "test_support.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <variant>

using namespace rl::bootstrap;
using namespace rl::bootstrap::test;
using rl::market::CurveInterpolation;
using Catch::Matchers::WithinAbs;

namespace {
/// |PV| of every instrument that has a PV (swaps, basis swaps, deposits, FRAs), per unit notional.
double maxAbsPv(const BootstrapProblem& p, const rl::math::Vector& z) {
    auto curves = p.curves(z);
    double worst = 0.0;
    for (const auto& bi : p.instruments()) {
        double pv = 0.0;
        if (auto* s = std::get_if<rl::instruments::Swap>(&bi.trade)) pv = rl::linear::pv(*s, curves, p.fixings());
        else if (auto* b = std::get_if<rl::instruments::BasisSwap>(&bi.trade)) pv = rl::linear::pv(*b, curves, p.fixings());
        else if (auto* d = std::get_if<rl::instruments::Deposit>(&bi.trade)) pv = rl::linear::pv(*d, curves);
        else if (auto* f = std::get_if<rl::instruments::Fra>(&bi.trade)) pv = rl::linear::pv(*f, curves);
        worst = std::max(worst, std::fabs(pv));
    }
    return worst;
}
} // namespace

TEST_CASE("Global solve reprices every instrument to 1e-10 of notional, under both interpolations", "[bootstrap][solve]") {
    for (auto interp : {CurveInterpolation::LinearLogDf, CurveInterpolation::MonotoneConvex}) {
        BootstrapProblem p = usdProblem(interp);
        REQUIRE(p.size() == p.instruments().size());
        BootstrapResult r = solve(p);
        INFO("interpolation " << rl::market::toString(interp) << ", iterations " << r.iterations);
        REQUIRE(r.converged);
        CHECK(r.iterations <= 15);
        CHECK(r.maxAbsBp < 1e-8);
        CHECK(maxAbsPv(p, r.z) < 1e-10);
        RepricingReport rep = reprice(p, r.z);
        CHECK(rep.maxAbsBp < 1e-8);
        CHECK(rep.rows.size() == p.size());
        for (const auto& row : rep.rows) {
            INFO(row.id);
            CHECK(std::fabs(row.residualBp) < 1e-8);
            // The model value in quoted units agrees with the quote to display precision.
            CHECK_THAT(row.model, WithinAbs(row.quote, 1e-9));
        }
        // Every future's implied price equals its quote: the convexity adjustment is applied.
        auto curves = p.curves(r.z);
        for (const auto& bi : p.instruments())
            if (auto* f = std::get_if<rl::instruments::Future>(&bi.trade))
                CHECK_THAT(rl::linear::impliedFuturePrice(*f, curves), WithinAbs(bi.quote, 1e-10));
    }
}

TEST_CASE("Newton converges quadratically from the flat guess and from a bad guess", "[bootstrap][solve]") {
    BootstrapProblem p = usdProblem();
    BootstrapOptions opt;
    opt.verbose = true;
    BootstrapResult r = solve(p, opt);
    REQUIRE(r.converged);
    CHECK(r.log.size() >= static_cast<std::size_t>(r.iterations));
    // From y = 0 (every df = 1) the line search must still get there.
    rl::math::Vector zero = rl::math::Vector::Zero(static_cast<Eigen::Index>(p.size()));
    CHECK(scaledResiduals(p, zero).cwiseAbs().maxCoeff() > 100.0); // wildly off
    CHECK(scaledResiduals(p, r.z).cwiseAbs().maxCoeff() < 1e-8);
    CHECK(scaledResiduals(p, initialGuess(p)).cwiseAbs().maxCoeff() < 200.0); // flat guess is within a couple of percent
}

TEST_CASE("Scaled residuals are commensurable: one bp of quote error reads as one bp everywhere", "[bootstrap][scaling]") {
    BootstrapProblem p = usdProblem();
    BootstrapResult r = solve(p);
    REQUIRE(r.converged);
    const rl::math::Vector q = p.quoteValues();
    for (Eigen::Index i = 0; i < q.size(); ++i) {
        const auto& bi = p.instruments()[static_cast<std::size_t>(i)];
        rl::math::Vector qb = q;
        // One basis point of the quote (a futures price moves by 0.01 per bp).
        qb[i] += bi.source.unit == rl::market::QuoteUnit::Price ? -0.01 : 1e-4;
        BootstrapProblem pb = p.withQuotes(qb);
        rl::math::Vector s = scaledResiduals(pb, r.z);
        INFO(bi.id);
        CHECK_THAT(std::fabs(s[i]), WithinAbs(1.0, 0.05)); // within 5% of one bp (rate-form: exact; PV-form: annuity-scaled)
    }
}

TEST_CASE("Sequential bootstrap leaves the cross-assigned basis swaps mispriced; the global solve does not", "[bootstrap][sequential]") {
    BootstrapProblem p = usdProblem();
    BootstrapResult seq = solveSequential(p);
    BootstrapResult glob = solve(p);
    REQUIRE(glob.converged);
    const rl::math::Vector sSeq = scaledResiduals(p, seq.z);
    const rl::math::Vector sGlob = scaledResiduals(p, glob.z);
    double worstSeq = 0.0;
    std::string worstId;
    for (std::size_t i = 0; i < p.size(); ++i) {
        const auto& bi = p.instruments()[i];
        const double e = std::fabs(sSeq[static_cast<Eigen::Index>(i)]);
        if (e > worstSeq) {
            worstSeq = e;
            worstId = bi.id;
        }
        // SOFR-only instruments (OIS on SOFR, futures) reprice in both.
        if (bi.curveName == "USD-SOFR" && bi.type != rl::market::QuoteType::BasisSwap) CHECK(e < 1e-6);
    }
    INFO("worst sequential residual " << worstSeq << " bp on " << worstId);
    CHECK(worstSeq > 0.1);
    CHECK(worstId.rfind("FF-SOFR-", 0) == 0);
    CHECK(sGlob.cwiseAbs().maxCoeff() < 1e-8);
}

TEST_CASE("Determinism and the snapshot hash", "[bootstrap][solve][hash]") {
    BootstrapProblem p = usdProblem();
    BootstrapResult a = solve(p);
    BootstrapResult b = solve(p);
    REQUIRE(a.converged);
    CHECK(a.z == b.z); // bit-identical
    auto s1 = rl::market::SnapshotBuilder<double>().asOf(p.asOf()).curves(p.curves(a.z)).fixings(p.fixings()).quotes(p.quotes()).build();
    auto s2 = rl::market::SnapshotBuilder<double>().asOf(p.asOf()).curves(p.curves(b.z)).fixings(p.fixings()).quotes(p.quotes()).build();
    CHECK(s1->hash() == s2->hash());
    BootstrapProblem pb = p.withQuotes(p.quoteValues() + rl::math::Vector::Constant(static_cast<Eigen::Index>(p.size()), 1e-6));
    BootstrapResult c = solve(pb);
    auto s3 = rl::market::SnapshotBuilder<double>().asOf(p.asOf()).curves(pb.curves(c.z)).fixings(p.fixings()).quotes(pb.quotes()).build();
    CHECK(s3->hash() != s1->hash());
}

TEST_CASE("Turn of year: the jump is two days wide", "[bootstrap][turn]") {
    // The USD file carries one 8 bp turn over 2026-12-31 → 2027-01-04; it becomes the last parameter of
    // the SOFR block (YieldCurve jump layer). Everything else on the curve must be unaffected.
    BootstrapProblem p = usdProblem(CurveInterpolation::LinearLogDf, true);
    REQUIRE(p.parameterCount() == p.instruments().size());
    REQUIRE(p.curveSpecs()[0].turns.size() == 1);
    BootstrapResult r = solve(p);
    REQUIRE(r.converged);
    auto curves = p.curves(r.z);
    const auto& sofr = curves.discount();
    auto fwd = [&](rl::core::Date a, rl::core::Date b) { return rl::linear::continuousForward(sofr, sofr.time(a), sofr.time(b)); };
    const double turn = fwd(rl::core::Date(2026, 12, 31), rl::core::Date(2027, 1, 4));
    const double before = fwd(rl::core::Date(2026, 12, 21), rl::core::Date(2026, 12, 28));
    const double after = fwd(rl::core::Date(2027, 1, 6), rl::core::Date(2027, 1, 13));
    CHECK_THAT(turn - 0.5 * (before + after), WithinAbs(8e-4, 1e-6));
    CHECK_THAT(before, WithinAbs(after, 5e-4)); // no spike leaking into the neighbours
}
