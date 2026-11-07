#pragma once
// The exposure of a Bermudan swaption held to its counterparty.
//
// A swap's exposure needs only today's model: its future value is a known function of the state.
// A Bermudan's does not have one: its future value is the value of an option that may or may not
// have been exercised already, and after exercise the holder owns a swap instead. On each path, at
// each date, the position is one of two things:
//
//     not yet exercised   the option, worth max(exercise value, continuation value)  ≥ 0
//     exercised earlier   the remaining swap, worth whatever rates have made it      either sign
//
// So the profile needs the continuation value as a FUNCTION of the state at every exercise date,
// which is the same object the pricer builds. Here it is computed by backward induction on a grid:
// between exercise dates the state moves by an exact Gaussian step, so
//
//     C_i(x) = ∫ V_{i+1}(x') φ(x' − x; ζ_{i+1} − ζ_i) dx',      V_i = max(h_i, C_i),
//
// with no time discretisation. Paths are then simulated forward and exercised when h_i ≥ C_i, the
// optimal rule, so the same simulation gives the price (checked against the PDE), the exercise
// probabilities, and the exposure.
//
// What to expect. Before the first exercise date the position is an option: exposure is its price,
// never negative. As paths exercise, option exposure is replaced by swap exposure, which is larger
// for the paths that exercised (they did so because the swap was valuable) but then amortises, and
// can go negative. Physical settlement is assumed; a cash-settled Bermudan has no exposure after
// exercise at all.
//
// The profile is sampled on the exercise dates, which for a co-terminal Bermudan are the coupon
// dates of the underlying swap; after exercise the remaining swap at date j is exercise right j.

#include "rl/pricers/bermudan.hpp"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace rl::xva {

struct BermudanExposureOptions {
    std::size_t paths = 50000;
    std::uint64_t seed = 20260904;
    int gridPoints = 801;    ///< state grid at each exercise date (odd, so that x = 0 is a node)
    double stdDevs = 8.0;
};

struct BermudanExposure {
    std::vector<double> times;            ///< 0, then the exercise dates
    std::vector<double> epe;              ///< discounted E[max(position, 0)]
    std::vector<double> ene;              ///< discounted E[min(position, 0)]: only from exercised paths
    std::vector<double> epeStandardError;
    std::vector<double> optionEpe;        ///< the part of EPE from paths still holding the option
    std::vector<double> exercisedBy;      ///< fraction of paths exercised on or before each date
    double price = 0.0;                   ///< today's value from the grid
    double simulatedPrice = 0.0;          ///< E[exercise value at the exercise date] over the paths
    double simulatedPriceStandardError = 0.0;
};

/// Exposure of a LONG position in the Bermudan (the holder's exposure to the writer).
BermudanExposure bermudanExposure(const pricers::BermudanProblem& problem, const BermudanExposureOptions& options = {});

} // namespace rl::xva
