#pragma once
// Monte Carlo exposure of a netting set of swaps, with and without collateral.
//
// Counterparty risk is about the value of a PORTFOLIO at future dates: if the counterparty defaults
// at t, what is lost is the positive part of the netted value of every trade under the same
// agreement, less whatever collateral is held. So the engine simulates the model's state forward,
// revalues every trade on every path at every exposure date, and aggregates:
//
//     EPE(t)   E[ max(V(t), 0) ]     expected positive exposure: what the counterparty owes us
//     ENE(t)   E[ min(V(t), 0) ]     expected negative exposure: what we owe (for DVA and funding)
//     PFE(t)   a high quantile of max(V(t), 0)
//     gross    Σ_trades E[ max(V_i(t), 0) ]: the same trades with no netting agreement
//
// Everything is in DEFLATED units, V/N, as in the pricers: the expectation of a deflated value is
// its price today, so EPE(t) is already discounted and CVA is EPE weighted by default probabilities.
//
// Under LGM the simulation is exact (the state at each date is the previous state plus an
// independent Gaussian increment) and a swap is a sum of bonds, each one exponential in the state.
//
// NETTING. V(t) = Σ_i V_i(t) inside the max. max(Σ V_i, 0) ≤ Σ max(V_i, 0), so netting can only
// reduce exposure; by how much depends on how the trades offset, which the engine measures.
//
// COLLATERAL. Under a credit support annex the party that is out of the money posts collateral
// against the portfolio's value above a THRESHOLD, in amounts of at least the MINIMUM TRANSFER
// AMOUNT. Collateral does not remove exposure, for one reason: time. After a default the surviving
// party needs the MARGIN PERIOD OF RISK (typically 10 business days) to notice, close out and
// re-hedge, and during it the portfolio moves while no more collateral arrives. So
//
//     collateral held at t = what the agreement called for on the value at t − MPR
//     exposure(t)         = max( V(t) − collateral, 0 )
//
// and even with a zero threshold what is left is the move of the portfolio over the margin period.
// That move includes cashflows: when a coupon rolls off inside the margin period the portfolio's
// value jumps by the net coupon while the collateral held still reflects the value before it, so
// collateralised exposure spikes on coupon dates.
//
// SIMPLIFICATIONS, each chosen so that the engine can be checked against a closed form:
//   * A floating coupon in progress is valued as if its rate had been fixed at the period start at
//     the then-prevailing term rate. Exact on coupon dates (where the closed-form swaption strip
//     applies) and off by the convexity of realised compounding in between.
//   * A coupon leaves the valuation at its accrual end, not on its payment date two days later.
//   * Collateral is cash in the deflated unit: its deflated value is held constant over the margin
//     period. Threshold and minimum transfer amount are symmetric. No initial margin, no haircuts.
//   * Single currency, swaps only.

#include "rl/instruments/swap.hpp"
#include "rl/market/curve_set.hpp"
#include "rl/models/g2pp.hpp"
#include "rl/models/lgm.hpp"

#include <cstddef>
#include <cstdint>
#include <optional>
#include <vector>

namespace rl::xva {

/// A credit support annex. Amounts are in the trades' notional units.
struct Csa {
    double threshold = 0.0;             ///< uncollateralised amount each way
    double minimumTransfer = 0.0;       ///< no call below this
    double marginPeriodYears = 10.0 / 252.0;
};

struct ExposureOptions {
    std::size_t paths = 20000;
    std::uint64_t seed = 20260904;
    double quantile = 0.975;
    bool keepPaths = false; ///< also return the netted value and the state on every path (for wrong-way risk)
};

struct ExposureProfile {
    std::vector<double> times;            ///< exposure dates in years; times[0] = 0
    std::vector<double> epe;              ///< netted, uncollateralised, discounted
    std::vector<double> ene;              ///< ≤ 0
    std::vector<double> epeStandardError;
    std::vector<double> pfe;              ///< `quantile` of the netted discounted exposure
    std::vector<double> grossEpe;         ///< sum of the trades' own EPEs: no netting agreement
    std::vector<double> collateralisedEpe; ///< empty when no CSA was given
    std::vector<double> collateralisedEne;
    double valueToday = 0.0;              ///< netted value of the trades on today's curve

    /// Only with `keepPaths`: per exposure date, per path, the netted deflated value and the model's
    /// first state variable. Entry 0 (today) holds the deterministic value and a zero state.
    std::vector<std::vector<double>> pathValue, pathState;
};

/// Exposure dates for a set of trades: today, every `stepYears` out to the last accrual end, and
/// every fixed-leg period start (where exposure steps as coupons roll off).
std::vector<double> exposureTimes(const std::vector<instruments::Swap>& trades, const market::CurveSet<double>& curves, double stepYears);

/// The fixed-leg period starts of one swap after today, preceded by 0: the dates on which the
/// closed-form swaption-strip exposure is defined.
std::vector<double> periodStartTimes(const instruments::Swap& swap, const market::CurveSet<double>& curves);

/// Simulate. `times` must start at 0 and be strictly increasing. Every swap must be an overnight
/// swap projected off the discount curve with no spread. Throws InputError otherwise.
ExposureProfile simulateExposure(const models::Lgm1F& model, const std::vector<instruments::Swap>& trades, const market::CurveSet<double>& curves,
                                 const std::vector<double>& times, const std::optional<Csa>& csa = std::nullopt, const ExposureOptions& options = {});

/// The same under the two-factor model. With two factors the trades of a netting set no longer move
/// together: swaps of different maturities load differently on the fast and the slow factor, so how
/// much they offset depends on the factor correlation, which no single-maturity swaption strip sees.
ExposureProfile simulateExposure(const models::G2pp& model, const std::vector<instruments::Swap>& trades, const market::CurveSet<double>& curves,
                                 const std::vector<double>& times, const std::optional<Csa>& csa = std::nullopt, const ExposureOptions& options = {});

} // namespace rl::xva
