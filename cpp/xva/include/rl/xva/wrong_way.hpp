#pragma once
// Wrong-way risk: when the counterparty is more likely to default exactly when it owes the most.
//
// CVA as computed in xva.hpp multiplies expected exposure by default probability, which is right
// only if the two are independent. They often are not. A counterparty whose credit deteriorates
// when rates rise, facing a bank that receives floating, defaults disproportionately in the states
// where the bank's swap is deep in the money: WRONG-way risk. The opposite dependence is right-way
// risk, and it lowers CVA. The exposure that matters is not the average exposure but the exposure
// CONDITIONAL ON DEFAULT.
//
// The model (in the spirit of Hull & White, "CVA and wrong-way risk", 2012). The hazard rate on
// each interval is a function of the rates model's state at the end of it,
//
//     λ_j(path) = a_j · exp( β z_j − ½ β² ),      z_j = x(t_j) / √ζ(t_j)  ~ N(0, 1),
//
// so a one-standard-deviation rise in the rates state multiplies the hazard by e^β. β > 0 means
// credit worsens as rates rise; β = 0 is independence. Path survival is
// S_j = exp(−Σ_{k≤j} λ_k Δt_k), and
//
//     CVA = (1 − R) · E[ Σ_j ½ (V⁺_{j−1} + V⁺_j) · (S_{j−1} − S_j) ].
//
// CALIBRATION. β is the one new parameter; the levels a_j are not free. Whatever the dependence,
// the model must still reprice the counterparty's CDS curve, so for each interval a_j is solved so
// that survival to t_j, as the market prices it, is unchanged:
//
//     E[ w_j · S_j ] = Q_market(t_j),      w_j = exp( −H(t_j) x_j − ½ H(t_j)² ζ(t_j) ).
//
// The weight w_j is the deflated bond maturing at t_j divided by its price: paths are simulated
// under the pricing model's numeraire, and the market's survival probability belongs to the
// t_j-forward measure. Without it the calibration would be to the wrong probabilities as soon as
// β ≠ 0. With β = 0 every a_j is the curve's own hazard rate and the result is the independent CVA.
//
// What this is not. The link is deterministic (one factor drives both), β is an input, not
// estimated, and default does not itself move rates. It measures how much CVA moves for a stated
// dependence; choosing the dependence is the hard part and is not attempted here.

#include "rl/models/lgm.hpp"
#include "rl/xva/credit.hpp"
#include "rl/xva/exposure_engine.hpp"

#include <vector>

namespace rl::xva {

struct WrongWayResult {
    double cva = 0.0;              ///< with the dependence
    double independentCva = 0.0;   ///< the same paths with β = 0
    double standardError = 0.0;    ///< of `cva`
    std::vector<double> times;
    std::vector<double> epe;                 ///< unconditional, as in the profile
    std::vector<double> epeGivenDefault;     ///< E[V⁺ | default in the interval ending here] (entry 0 unused)
    std::vector<double> hazardLevel;         ///< a_j per interval
    std::vector<double> survivalError;       ///< E[w_j S_j] − Q(t_j): zero to the root-finder's tolerance
};

/// `profile` must have been simulated under `model` with `keepPaths` and no collateral.
/// Throws InputError otherwise, and CalibrationError if the survival curve cannot be matched.
WrongWayResult wrongWayCva(const ExposureProfile& profile, const models::Lgm1F& model, const HazardCurve& hazard, double recovery, double beta);

} // namespace rl::xva
