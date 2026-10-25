// The two Monte Carlo bounds against the PDE: lower <= PDE <= upper, within sampling error.
#include "rl/pricers/bermudan_mc.hpp"

#include "rl/pricers/bermudan_pde.hpp"
#include "test_support.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <algorithm>
#include <cmath>

using namespace rl::pricers;
using namespace rl::pricers::test;
using Catch::Matchers::WithinAbs;

namespace {
double pdeReference(const BermudanProblem& p) {
    PdeOptions o;
    o.spaceSteps = 1600;
    o.timeSteps = 1600;
    return priceBermudanPde(p, o).value;
}
} // namespace

TEST_CASE("LSMC lower bound <= PDE <= Andersen-Broadie upper bound, and the gap is small", "[pricers][mc][bounds]") {
    for (bool payer : {true, false}) {
        const CalibratedBermudan b = bermudan10nc1(0.035, payer);
        const double pde = pdeReference(b.problem);
        const LsmcResult lower = priceBermudanLsmc(b.problem);
        const UpperBoundResult upper = priceBermudanUpperBound(b.problem, lower);
        INFO("lower " << lower.lower << " +/- " << lower.standardError << "  pde " << pde << "  upper " << upper.upper << " +/- " << upper.standardError);
        CHECK(lower.lower <= pde + 3.0 * lower.standardError);
        CHECK(upper.upper >= pde - 3.0 * upper.standardError);
        // A cubic in the state is nearly the optimal rule in a one-factor model: both bounds within 1% of the PDE.
        CHECK(pde - lower.lower < 0.01 * pde);
        CHECK(upper.upper - pde < 0.01 * pde);
        // Exercise probabilities are probabilities.
        double total = 0.0;
        for (double q : lower.exerciseProbability) {
            CHECK(q >= 0.0);
            total += q;
        }
        CHECK(total <= 1.0);
        CHECK(total > 0.3);
    }
}

TEST_CASE("In-sample LSMC is biased high relative to the out-of-sample estimate", "[pricers][mc][bias]") {
    // With few training paths and a flexible basis the rule overfits its own paths.
    const CalibratedBermudan b = bermudan10nc1();
    double inSample = 0.0, outOfSample = 0.0;
    const int trials = 20;
    for (int k = 0; k < trials; ++k) {
        LsmcOptions o;
        o.trainingPaths = 500;
        o.pricingPaths = 20000;
        o.degree = 6;
        o.seed = 1000 + static_cast<std::uint64_t>(k);
        const LsmcResult r = priceBermudanLsmc(b.problem, o);
        inSample += r.inSample / trials;
        outOfSample += r.lower / trials;
    }
    INFO("in-sample " << inSample << " out-of-sample " << outOfSample);
    CHECK(inSample > outOfSample);
    CHECK(outOfSample < pdeReference(b.problem));
}

TEST_CASE("A poor exercise rule gives a lower lower bound and a wider duality gap", "[pricers][mc][gap]") {
    const CalibratedBermudan b = bermudan10nc1();
    LsmcOptions good, poor;
    poor.degree = 0; // continuation value = a constant: the rule cannot see the state
    const LsmcResult lg = priceBermudanLsmc(b.problem, good), lp = priceBermudanLsmc(b.problem, poor);
    const UpperBoundResult ug = priceBermudanUpperBound(b.problem, lg), up = priceBermudanUpperBound(b.problem, lp);
    const double pde = pdeReference(b.problem);
    CHECK(lp.lower < lg.lower);
    CHECK(up.gap > 2.0 * ug.gap);
    CHECK(ug.gap > -3.0 * ug.gapStandardError); // the gap is non-negative up to sampling noise
    // Both are still valid brackets of the true price.
    CHECK(lp.lower <= pde + 3.0 * lp.standardError);
    CHECK(up.upper >= pde - 3.0 * up.standardError);
}

TEST_CASE("Single exercise date: LSMC is plain Monte Carlo of the European", "[pricers][mc]") {
    const CalibratedBermudan b = bermudan10nc1();
    const BermudanProblem one = subset(b.problem, {4});
    const LsmcResult r = priceBermudanLsmc(one);
    CHECK_THAT(r.lower, WithinAbs(one.europeanValues()[0], 4.0 * r.standardError));
    // Reproducible: same seed, same number.
    CHECK(priceBermudanLsmc(one).lower == r.lower);
}

TEST_CASE("Both regression bases give valid, close lower bounds", "[pricers][mc][basis]") {
    const CalibratedBermudan b = bermudan10nc1();
    LsmcOptions s, v;
    v.basis = RegressionBasis::ExerciseValue;
    const LsmcResult rs = priceBermudanLsmc(b.problem, s), rv = priceBermudanLsmc(b.problem, v);
    const double pde = pdeReference(b.problem);
    CHECK(rs.lower <= pde + 3.0 * rs.standardError);
    CHECK(rv.lower <= pde + 3.0 * rv.standardError);
    CHECK(std::fabs(rs.lower - rv.lower) < 0.005 * pde);
}
