#pragma once
// Bermudan swaption by Monte Carlo: a lower bound and an upper bound.
//
// Monte Carlo cannot price an early-exercise option directly, because the exercise decision needs
// the continuation value, which is itself an expectation. Two estimators bracket the price.
//
// LOWER BOUND — Longstaff & Schwartz (2001). Estimate the continuation value at each date by
// regressing the realised future payoff on functions of the state, across paths. That defines an
// exercise rule: exercise when the exercise value beats the fitted continuation. ANY exercise rule
// is a strategy the holder could follow, so its expected payoff is at most the true price. Two
// details decide whether the number means anything:
//   * the rule must be evaluated on paths INDEPENDENT of the ones it was fitted on. On its own
//     training paths the regression has seen the future of each path and the estimate is biased
//     high; `inSample` reports that number so the bias is visible.
//   * the regression uses only in-the-money paths: out of the money the decision is "do not
//     exercise" whatever the continuation value, so those paths only dilute the fit where it matters.
//
// UPPER BOUND — Andersen & Broadie (2004), from the dual formulation of Rogers / Haugh–Kogan:
//     price ≤ E[ max_i ( h_i − M_i ) ]     for any martingale M with M_0 = 0,
// with equality for the martingale part of the true value process. Take M from the lower bound's
// exercise rule: its increments are L_i − E_{i−1}[L_i], where L is the value of following the rule,
// and the conditional expectations are estimated by nested simulation. The better the rule, the
// tighter the bound; the difference upper − lower (the duality gap) is a direct measure of how much
// value the exercise rule leaves on the table.
//
// Under LGM both simulations are exact: the state at exercise date i is the previous state plus an
// independent N(0, ζ_i − ζ_{i−1}) draw, with no time-stepping error, and deflated payoffs need no
// discounting.

#include "rl/pricers/bermudan.hpp"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace rl::pricers {

enum class RegressionBasis {
    State,         ///< powers of the standardised state x/√ζ
    ExerciseValue, ///< powers of the exercise value itself (the swap's value)
};

/// An exercise rule: at date i, exercise iff the exercise value is positive and at least the fitted
/// continuation value. At the last date, exercise iff positive.
struct ExercisePolicy {
    RegressionBasis basis = RegressionBasis::State;
    int degree = 3;
    std::vector<std::vector<double>> coefficients; ///< per exercise date; empty = "do not exercise here"
    std::vector<double> scale;                     ///< per date: normalisation of the regressor

    double continuation(const BermudanProblem& problem, std::size_t i, double x, double exerciseValue) const;
    bool exercise(const BermudanProblem& problem, std::size_t i, double x, double exerciseValue) const;
};

struct LsmcOptions {
    std::size_t trainingPaths = 20000;
    std::size_t pricingPaths = 100000; ///< independent of the training paths
    RegressionBasis basis = RegressionBasis::State;
    int degree = 3;
    bool antithetic = true;            ///< pricing pass only
    std::uint64_t seed = 20260904;
};

struct LsmcResult {
    double lower = 0.0;         ///< out-of-sample estimate: a valid lower bound (up to its standard error)
    double standardError = 0.0;
    double inSample = 0.0;      ///< the same rule on its own training paths: biased high
    std::vector<double> exerciseProbability; ///< fraction of pricing paths exercised at each date
    ExercisePolicy policy;
};

LsmcResult priceBermudanLsmc(const BermudanProblem& problem, const LsmcOptions& options = {});

struct UpperBoundOptions {
    std::size_t outerPaths = 2000;
    std::size_t innerPaths = 500; ///< per outer path per exercise date
    std::uint64_t seed = 77;
};

struct UpperBoundResult {
    double gap = 0.0;              ///< Δ = E[max_i (h_i − M_i)] with M_0 = lower bound: the duality gap, ≥ 0 up to noise
    double gapStandardError = 0.0;
    double upper = 0.0;            ///< lower + gap
    double standardError = 0.0;    ///< of `upper`: the lower bound's and the gap's errors combined
};

/// Andersen–Broadie upper bound from the exercise rule in `lower`. The martingale starts at the
/// rule's value today, so the upper bound is lower + gap; the gap itself is estimated far more
/// precisely than either bound, because along each path max_i (h_i − M_i) is nearly constant when
/// the rule is nearly optimal.
UpperBoundResult priceBermudanUpperBound(const BermudanProblem& problem, const LsmcResult& lower, const UpperBoundOptions& options = {});

} // namespace rl::pricers
