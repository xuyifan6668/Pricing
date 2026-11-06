#pragma once
// From exposure to price adjustments.
//
// A derivative's model price assumes both parties perform and that cash can be borrowed and lent
// at the risk-free rate. The adjustments put those assumptions back:
//
//   CVA   the expected loss if the COUNTERPARTY defaults while owing us:
//             CVA = (1 − R_c) · Σ_i EPE(t_i) · P(counterparty defaults in (t_{i−1}, t_i])
//   DVA   the mirror image: the expected gain to us if WE default while owing them:
//             DVA = (1 − R_own) · Σ_i |ENE(t_i)| · P(we default in (t_{i−1}, t_i])
//   FVA   the cost of funding. An uncollateralised asset has to be funded at our borrowing spread
//         over the risk-free rate; an uncollateralised liability provides funding:
//             FCA = s_f · ∫ EPE(t) dt,     FBA = s_f · ∫ |ENE(t)| dt,     FVA = FCA − FBA.
//
// EPE and ENE are discounted expectations (exposure_engine.hpp), so no further discounting appears.
// Each sum uses the average of the exposure at the two ends of its interval.
//
// SIMPLIFICATIONS. Default and exposure are independent: no wrong-way risk. CVA does not condition
// on our own survival, nor DVA on the counterparty's (unilateral adjustments). FVA is not weighted
// by survival, and FBA overlaps with DVA, as it does in practice: the two are reported separately
// and not added. A single flat funding spread.
//
// For one swap under LGM no simulation is needed: EPE(t) at a coupon date is the price of the
// European swaption on the remaining swap, in closed form. That gives a CVA with no Monte Carlo
// error against which the engine is checked.

#include "rl/instruments/swap.hpp"
#include "rl/market/curve_set.hpp"
#include "rl/models/lgm.hpp"
#include "rl/xva/credit.hpp"

#include <vector>

namespace rl::xva {

struct XvaInputs {
    HazardCurve counterparty = HazardCurve::flat(0.0);
    double counterpartyRecovery = 0.4;
    HazardCurve own = HazardCurve::flat(0.0);
    double ownRecovery = 0.4;
    double fundingSpread = 0.0; ///< decimal, flat
};

struct XvaResult {
    double cva = 0.0; ///< ≥ 0: a cost
    double dva = 0.0; ///< ≥ 0: a benefit
    double fca = 0.0; ///< ≥ 0: a cost
    double fba = 0.0; ///< ≥ 0: a benefit
    double fva() const noexcept { return fca - fba; }
};

/// (1 − recovery) · Σ ½ (e_{i−1} + e_i) · P(default in (t_{i−1}, t_i]) for a non-negative exposure profile.
double creditAdjustment(const std::vector<double>& times, const std::vector<double>& exposure, const HazardCurve& hazard, double recovery);

/// spread · Σ ½ (e_{i−1} + e_i) (t_i − t_{i−1}).
double fundingAdjustment(const std::vector<double>& times, const std::vector<double>& exposure, double spread);

/// All four adjustments from a profile's EPE and ENE (ENE ≤ 0).
XvaResult computeXva(const std::vector<double>& times, const std::vector<double>& epe, const std::vector<double>& ene, const XvaInputs& inputs);

/// EPE and ENE of one swap at its fixed-leg period starts, in closed form: a strip of European
/// swaptions on the remaining swap (Sorensen & Bollier). times[0] = 0.
struct SwapExposure {
    std::vector<double> times, epe, ene;
};
SwapExposure swapExposureClosedForm(const models::Lgm1F& model, const instruments::Swap& swap, const market::CurveSet<double>& curves);

} // namespace rl::xva
