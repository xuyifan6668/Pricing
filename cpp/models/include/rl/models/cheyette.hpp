#pragma once
// The one-factor quasi-Gaussian (Cheyette) model with local volatility: a short-rate model that has
// a smile.
//
// Why. Hull-White / LGM has one volatility per date, whatever the level of rates. Its swaption smile
// is therefore almost flat in normal vol, and it can be made to agree with the market at one strike
// per expiry only (pricers/bermudan.hpp, calibrateBermudanAtStrikes). The market's vol rises with
// the strike here, and a Bermudan is exercised well away from its own strike. Letting the short
// rate's volatility depend on the rate itself gives the model a skew to calibrate.
//
// The model (Cheyette 1992; Ritchken & Sankarasubramanian 1995; Andersen & Piterbarg ch. 13). Under
// the risk-neutral measure, with r(t) = f(0, t) + x(t):
//
//     dx = (y − κ x) dt + η(t, x) dW,        x(0) = 0
//     dy = (η(t, x)² − 2 κ y) dt,            y(0) = 0
//
// y is the accumulated variance of x, discounted by mean reversion. It is not a second source of
// risk (there is one Brownian motion) but it is a second STATE: with a volatility that depends on
// the path, the variance accumulated so far is no longer a function of time alone, and the model
// needs to carry it to stay arbitrage-free. That is the price of the smile. The zero-coupon bond is
// still closed form in the two states:
//
//     P(t, T) = P(0, T)/P(0, t) · exp( −G(T − t) x − ½ G(T − t)² y ),      G(τ) = (1 − e^{−κτ})/κ.
//
// With η independent of x the variance y is deterministic and this is exactly Hull-White with
// short-rate volatility η(t): the tests check prices against the LGM closed form in that case.
//
// Local volatility. η(t, x) = λ(t) · (1 + b x), floored and capped. λ(t) is piecewise constant and
// sets the LEVEL of volatility through time; b sets the SKEW: b > 0 makes volatility rise with
// rates (toward lognormal behaviour), b < 0 fall. One number for the skew is the simplest choice
// that can be calibrated; a desk would let b depend on time and add stochastic volatility for the
// curvature of the smile.
//
// There are no closed forms for swaptions, so prices come from Monte Carlo
// (pricers/cheyette_mc.hpp). Immutable after construction.

#include <cstddef>
#include <vector>

namespace rl::models {

class Cheyette {
public:
    /// λ_i applies on (t_{i−1}, t_i], t_{−1} = 0, and the last continues beyond the last time.
    /// Throws InputError unless κ and the skew are finite, the times are positive and strictly
    /// increasing and the levels positive.
    Cheyette(double kappa, double skew, std::vector<double> times, std::vector<double> levels);

    double kappa() const noexcept { return kappa_; }
    double skew() const noexcept { return skew_; }
    const std::vector<double>& times() const noexcept { return times_; }
    const std::vector<double>& levels() const noexcept { return levels_; }

    /// λ(t): the level on the interval containing t (left-open, right-closed).
    double level(double t) const noexcept;
    /// η(t, x) = λ(t) · clamp(1 + b x, 0.05, 20).
    double localVolatility(double t, double x) const noexcept { return level(t) * shape(x); }
    /// The state-dependent factor alone.
    double shape(double x) const noexcept;
    /// G(τ) = (1 − e^{−κτ})/κ, and τ in the limit κ → 0.
    double G(double tau) const noexcept;
    /// P(t, T) · P(0, t)/P(0, T) = exp(−G(T − t) x − ½ G(T − t)² y): the bond relative to its forward price today.
    double bondRatio(double t, double T, double x, double y) const noexcept;

private:
    double kappa_, skew_;
    std::vector<double> times_, levels_;
};

} // namespace rl::models
