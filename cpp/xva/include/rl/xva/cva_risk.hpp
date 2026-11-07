#pragma once
// CVA sensitivities: what a CVA desk hedges with.
//
//   CS01      to each CDS spread the credit curve is built from: bump one quote by 1 bp, re-bootstrap
//             the hazard curve, recompute. Exposure does not change, so this costs nothing.
//
//   IR delta  to each par quote the discount curve is built from. Exposure depends on the curve, so
//             the naive route is one full CVA per quote. For one swap, whose CVA is a strip of
//             swaptions in closed form, there is an adjoint route:
//
//               1. Each swaption is Σ_k c_k Φ_k, where c_k = (cashflow amount) × (discount factor)
//                  and Φ_k is the probability-like weight of bond k. Because the exercise boundary
//                  is optimal, the price's sensitivity to the boundary is zero, so
//                  ∂(swaption)/∂c_k = Φ_k exactly: the weights the price was computed with are
//                  already its derivatives.
//               2. CVA is therefore linear in the c_k with known coefficients. Record that one
//                  linear combination on the AAD tape as a function of the curve nodes and sweep
//                  once: sensitivity to every node.
//               3. One transposed solve with the bootstrap Jacobian (implicit-function theorem)
//                  turns node risk into par-quote risk.
//
//             The adjoint route holds the model (H, ζ), the credit curve and the dates fixed. That
//             is a choice with a meaning: it is the delta of a desk that hedges vega separately. The
//             bump route can do the same, or re-calibrate the model to the swaption vols after each
//             bump, which moves ζ and gives a different number. Both are reported.

#include "rl/instruments/swap.hpp"
#include "rl/models/lgm.hpp"
#include "rl/models/usd_market.hpp"
#include "rl/pricers/risk.hpp"
#include "rl/xva/credit.hpp"

#include <functional>
#include <vector>

namespace rl::xva {

/// CVA of one swap from the closed-form swaption strip.
double swapCvaClosedForm(const models::Lgm1F& model, const instruments::Swap& swap, const market::CurveSet<double>& curves, const HazardCurve& hazard,
                         double recovery);

/// Change in a credit adjustment for +1 bp on each CDS quote in turn. `adjustment(curve)` computes
/// the adjustment on a given hazard curve.
std::vector<pricers::RiskRow> cs01Ladder(const std::vector<CdsQuote>& quotes, double recovery, const std::function<double(double)>& discount,
                                         const std::function<double(const HazardCurve&)>& adjustment, double bumpBp = 1.0);

/// Par-quote delta of one swap's closed-form CVA by the adjoint route, model and credit curve held fixed.
std::vector<pricers::RiskRow> swapCvaDeltaAdjoint(const models::UsdMarket& market, const models::Lgm1F& model, const instruments::Swap& swap,
                                                  const HazardCurve& hazard, double recovery);

} // namespace rl::xva
