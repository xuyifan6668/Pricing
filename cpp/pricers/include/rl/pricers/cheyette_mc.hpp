#pragma once
// Monte Carlo for the Cheyette model: European swaptions, calibration of the volatility level to a
// co-terminal strip, and the Bermudan by least-squares Monte Carlo.
//
// Simulation. Under the risk-neutral measure the numeraire is the bank account, so the value today
// of a payoff V at t is E[ V · exp(−∫₀ᵗ r) ], and with r = f(0, ·) + x the discount factor is
// P(0, t) · exp(−∫₀ᵗ x). A swap's value at t is a sum of bonds P(t, T_k), each P(0, T_k)/P(0, t)
// times a known function of (x, y), so P(0, t) cancels and the discounted exercise value is
//
//     exp(−∫₀ᵗ x) · Σ_k a_k D(T_k) exp( −G(T_k − t) x − ½ G(T_k − t)² y ).
//
// (x, y) are stepped by Euler with a fixed number of steps a year and ∫x by the trapezoid rule.
// Unlike the LGM pricers this has a time-discretisation error; with the skew switched off the model
// is Hull-White and the error is measured against the closed form (tests, RESULTS).
//
// Variance. Paths come in antithetic pairs, which removes most of the noise of a payoff that is
// close to linear in the driver. Comparisons between two skews use the same random numbers, so the
// DIFFERENCE between two prices is far more precise than either price.
//
// Calibration. The level λ on each interval is solved in turn so that the swaption expiring at its
// end is repriced, as ζ is in LGM. The price is a Monte Carlo estimate, so each interval's paths
// are re-simulated from the stored state at its start with the same random numbers for every trial
// λ: the objective is then a smooth function of λ and a root-finder applies. The calibration is
// exact on its own paths and accurate on fresh ones to their sampling error.

#include "rl/models/cheyette.hpp"
#include "rl/models/flows.hpp"
#include "rl/models/lgm_calibration.hpp"

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace rl::pricers {

struct CheyetteOptions {
    std::size_t paths = 100000; ///< even: antithetic pairs
    int stepsPerYear = 24;
    std::uint64_t seed = 20260904;
};

struct CheyettePrice {
    double value = 0.0;
    double standardError = 0.0;
};

/// European options to enter each of `rights` (any set, sorted by exercise time; several may share a
/// time, e.g. the same swaption at different strikes). One simulation prices them all.
std::vector<CheyettePrice> priceEuropeansCheyette(const models::Cheyette& model, const std::vector<models::ExerciseRight>& rights, const CheyetteOptions& options = {});

struct CheyetteCalibrationRow {
    std::string label;
    double expiry = 0.0;
    double marketVolBp = 0.0, modelVolBp = 0.0, errorBp = 0.0;
    double levelBp = 0.0; ///< λ on the interval ending at this expiry, bp
};

struct CheyetteCalibration {
    models::Cheyette model;
    std::vector<CheyetteCalibrationRow> rows;
    double maxAbsErrorBp = 0.0; ///< on the calibration's own paths
    bool ok = false;
};

/// Bootstrap λ to a strip in increasing expiry order, for a given mean reversion and skew. Throws
/// InputError on an empty or unordered strip and CalibrationError if a price is not attainable.
CheyetteCalibration calibrateCheyette(double kappa, double skew, const std::vector<models::CalibrationInstrument>& strip, const CheyetteOptions& options = {});

struct CheyetteLsmcOptions {
    std::size_t trainingPaths = 50000;
    CheyetteOptions pricing;  ///< pricing pass; its seed is not used for training
    int degree = 3;
};

struct CheyetteLsmcResult {
    double lower = 0.0; ///< out-of-sample: a lower bound up to its standard error and the time-stepping error
    double standardError = 0.0;
    std::vector<double> exerciseProbability;
};

/// The Bermudan with these exercise rights (increasing in time). Regression on the short-rate state
/// x and on the accumulated variance y, in-the-money paths only.
CheyetteLsmcResult priceBermudanCheyette(const models::Cheyette& model, const std::vector<models::ExerciseRight>& rights, const CheyetteLsmcOptions& options = {});

} // namespace rl::pricers
