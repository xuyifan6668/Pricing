#pragma once
// Risk: what moves the price, expressed in the instruments a desk hedges with.
//
//   Curve delta   sensitivity to each par quote the curve is built from (swap rates, futures prices,
//                 basis spreads). Two routes to the same numbers:
//                   bump     shift one quote by 1 bp, re-bootstrap, reprice; one bootstrap per quote.
//                   adjoint  one reverse sweep through the pricing code gives the sensitivity to every
//                            curve node at once; one transposed linear solve with the bootstrap
//                            Jacobian (implicit-function theorem) turns node risk into par-quote risk.
//                 The adjoint route is available for anything priced directly off the curve (swaps).
//   Model prices  (the Bermudan) are bumped: shift a quote, rebuild the market, RE-CALIBRATE the
//                 model, reprice. Re-calibrating matters: the hedge instruments are the co-terminal
//                 swaptions, and the model must keep repricing them after the bump.
//   Vega          sensitivity to each co-terminal swaption's normal vol (1 bp), by re-calibrating
//                 with that one instrument's vol bumped. This is the vega a Bermudan desk reads: it
//                 says how many of each European to trade against the Bermudan.

#include "rl/instruments/swap.hpp"
#include "rl/models/usd_market.hpp"
#include "rl/pricers/bermudan.hpp"

#include <functional>
#include <string>
#include <vector>

namespace rl::pricers {

struct RiskRow {
    std::string id;
    double value = 0.0; ///< price change for +1 bp of the quote (rate up 1 bp; a future's price DOWN 0.01)
};

/// Par-quote delta of a portfolio of swaps by the adjoint route. PVs are per the swaps' notionals.
std::vector<RiskRow> swapParDeltaAdjoint(const models::UsdMarket& market, const std::vector<instruments::Swap>& swaps);

/// The same by bump-and-rebuild: one bootstrap per quote. The quote is moved by `bumpBp` basis points
/// and the result scaled back to "per 1 bp".
std::vector<RiskRow> swapParDeltaBump(const models::UsdMarket& market, const std::vector<instruments::Swap>& swaps, double bumpBp = 1.0);

/// Par-quote delta of any price that is a function of the market, by bump-and-rebuild (+1 bp, one-sided).
std::vector<RiskRow> parDeltaBump(const models::UsdMarket& market, const std::function<double(const models::UsdMarket&)>& price, double bumpBp = 1.0);

/// Vega of a Bermudan to each co-terminal calibration swaption: bump that instrument's normal vol
/// by 1 bp, re-bootstrap ζ, reprice with `price`.
std::vector<RiskRow> bermudanVegaByExpiry(const models::UsdMarket& market, const instruments::Swap& swap, double firstCallYears, int noticeDays, double kappa,
                                          const std::function<double(const BermudanProblem&)>& price);

} // namespace rl::pricers
