#pragma once
// Counterparty exposure of a swap under the same model that prices the Bermudan.
//
// The expected positive exposure at a future date t is EE(t) = E[ max(V_swap(t), 0) ] (discounted to
// today). That is, by definition, the price of a European option to enter the remaining swap at t:
// a swap's exposure profile IS a strip of co-terminal swaptions (Sorensen & Bollier 1994), the same
// instruments the Bermudan is calibrated to. So the front-office model gives the exposure profile in
// closed form, and a Monte Carlo exposure engine can be validated against it, which is what the
// test does.
//
// The Bermudan connects to this directly: a counterparty holding the right to cancel a swap holds a
// Bermudan on the offsetting swap, and CVA on a swap is a credit-weighted sum of these swaptions.
//
// SIMPLIFICATION: exposure is sampled on coupon period starts only (no accrued interest between
// dates), with no collateral, netting or margin period of risk.

#include "rl/instruments/swap.hpp"
#include "rl/market/curve_set.hpp"
#include "rl/models/lgm.hpp"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace rl::pricers {

struct ExposurePoint {
    double time = 0.0;
    double expectedExposure = 0.0;   ///< discounted EE: closed form (European swaption on the remaining swap)
    double expectedExposureMc = 0.0; ///< the same by simulation
    double mcStandardError = 0.0;
    double potentialExposure = 0.0;  ///< discounted exposure at the `quantile` state (closed form: the swap is monotone in x)
    double expectedNegativeExposure = 0.0; ///< discounted E[min(V, 0)], closed form (≤ 0)
};

struct ExposureOptions {
    std::size_t paths = 50000;
    double quantile = 0.975;
    std::uint64_t seed = 11;
};

/// Exposure profile of `swap` (from its holder's side) at each fixed period start after today.
std::vector<ExposurePoint> swapExposureProfile(const models::Lgm1F& model, const instruments::Swap& swap, const market::CurveSet<double>& curves,
                                               const ExposureOptions& options = {});

} // namespace rl::pricers
