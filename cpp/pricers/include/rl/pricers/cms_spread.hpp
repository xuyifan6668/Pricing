#pragma once
// A CMS spread option under the two-factor model: the product that prices correlation.
//
// The option pays, at expiry T,   max( S_long(T) − S_short(T) − K, 0 )   (a call; a put pays the
// reverse), where S_long and S_short are the par rates then of two swaps of different tenors, say
// 10 years and 2 years. It is a bet on the SLOPE of the curve.
//
// Why it needs two factors. The spread's variance is
//     Var(S_long − S_short) = σ_long² + σ_short² − 2 ρ_rates σ_long σ_short.
// The two volatilities are swaption vols, which any model can be calibrated to. The correlation
// ρ_rates between the two rates is not in any swaption price. In a one-factor model it is very
// nearly 1 whatever the parameters, so once both vols are matched the spread has almost no variance
// left and the option is worth close to its intrinsic value. A two-factor model can set ρ_rates,
// and the option's price is a direct reading of it. That is why CMS spread options are what a desk
// marks the factor correlation to.
//
// Pricing. Under the measure whose numeraire is the bond maturing at T the state is Gaussian with
// no drift (g2pp.hpp), and each swap rate at T is a ratio of sums of bonds, each an exponential in
// the state. So
//     price = D(T) · E[ payoff(x) ],     x ~ N(0, C(T)),
// a two-dimensional integral of a kinked function, done by Simpson's rule in two independent
// standard normals. The same integral gives the means, volatilities and correlation of the two
// rates, which is what makes the price interpretable.
//
// SIMPLIFICATION: one payment at expiry with unit accrual, and no payment lag after the fixing. A
// traded CMS spread option is a strip of such caplets with accrual fractions; each is this payoff.
// The swap rates are the model's own par rates of the two swaps, convexity included.

#include "rl/models/flows.hpp"
#include "rl/models/g2pp.hpp"

namespace rl::pricers {

struct CmsSpreadResult {
    double price = 0.0;            ///< today's value per unit notional
    double forwardSpread = 0.0;    ///< E[S_long − S_short] under the T-forward measure
    double longRate = 0.0;         ///< E[S_long]: the convexity-adjusted forward CMS rate
    double shortRate = 0.0;
    double longVolBp = 0.0;        ///< standard deviation of S_long at expiry, annualised, bp
    double shortVolBp = 0.0;
    double spreadVolBp = 0.0;      ///< the same for the spread
    double rateCorrelation = 0.0;  ///< correlation between the two rates at expiry
};

/// `longSwap` and `shortSwap` are the rights to enter the two swaps at the same exercise time (their
/// strikes are irrelevant: only their cashflow dates are used). `dfExpiry` is D(T). `pointsPerAxis`
/// must be odd. Throws InputError if the two rights have different exercise times.
CmsSpreadResult priceCmsSpread(const models::G2pp& model, const models::ExerciseRight& longSwap, const models::ExerciseRight& shortSwap, double strike,
                               bool call, double dfExpiry, int pointsPerAxis = 401);

} // namespace rl::pricers
