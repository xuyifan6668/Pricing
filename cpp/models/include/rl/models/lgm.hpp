#pragma once
// The one-factor Linear Gaussian Markov model (Hagan, "Evaluating and hedging exotic swap
// instruments via LGM"). It is Hull–White in a parameterisation that makes pricing and calibration
// simpler; the two are the same model.
//
// One state variable x with x_0 = 0 and, under the model's numeraire measure,
//
//     dx_t = sqrt(ζ'(t)) dW_t          ⇒          x_t ~ N(0, ζ(t)),   no drift.
//
// The numeraire and the zero-coupon bond are
//
//     N(t, x)    = exp( H(t) x + ½ H(t)² ζ(t) ) / D(t)
//     P(t, T; x) = D(T)/D(t) · exp( −(H(T) − H(t)) x − ½ (H(T)² − H(t)²) ζ(t) )
//
// with D(t) = P(0, t) today's discount curve. So the DEFLATED bond is a single exponential:
//
//     P(t, T; x) / N(t, x) = D(T) · exp( −H(T) x − ½ H(T)² ζ(t) ).
//
// Three consequences that the pricers rely on:
//   1. The model fits today's curve by construction: D(t) is an input, not a calibrated quantity.
//   2. The deflated value of any payoff is a martingale in a driftless Gaussian state, so a backward
//      PDE in the "clock" ζ is the heat equation, and Monte Carlo steps are exact Gaussian draws.
//   3. A European swaption expiring at t depends on the model only through ζ(t) and H at its cashflow
//      dates. With H fixed, a strip of swaptions with increasing expiries pins ζ one point at a
//      time: calibration is a sequence of 1-D root searches, not a global optimisation.
//
// Parameters and their Hull–White meaning (short-rate dr = (θ(t) − κ r) dt + σ(t) dW):
//     H(t) = (1 − e^{−κt}) / κ                     κ is the mean reversion
//     ζ(t) = ∫₀ᵗ σ(s)² e^{2κs} ds                  σ(t) is the short-rate volatility
// σ(t) is piecewise constant between the ζ nodes, which is what the interpolation of ζ implements.
//
// What κ does. κ does not change the price of the calibration swaptions (ζ absorbs it). It changes
// how variance is distributed ACROSS time: with higher κ, a shock to rates today has decayed by the
// time later exercise dates arrive, so rates at different exercise dates are less correlated. A
// Bermudan is an option on the best of several correlated exercise decisions, so its value above the
// best European (the "switch value") rises with κ. κ is therefore the parameter a Bermudan desk
// actually argues about; here it is an input.
//
// Immutable after construction; safe to share across threads.

#include <cstddef>
#include <vector>

namespace rl::models {

class Lgm1F {
public:
    /// ζ nodes at strictly increasing times > 0 with non-decreasing values > 0; ζ(0) = 0 is implicit.
    /// Throws InputError otherwise, or if κ is not finite.
    Lgm1F(double kappa, std::vector<double> zetaTimes, std::vector<double> zetaValues);

    /// A model with constant Hull–White short-rate volatility σ out to `horizon` years.
    static Lgm1F constantSigma(double kappa, double sigma, double horizon);

    double kappa() const noexcept { return kappa_; }
    const std::vector<double>& zetaTimes() const noexcept { return times_; }
    const std::vector<double>& zetaValues() const noexcept { return zetas_; }

    /// H(t) = (1 − e^{−κt})/κ, and t in the limit κ → 0.
    double H(double t) const noexcept { return H(kappa_, t); }
    static double H(double kappa, double t) noexcept;

    /// ζ(t): the variance of the state at t. Interpolated so that the Hull–White σ is constant on
    /// each interval; beyond the last node the last interval's σ continues.
    double zeta(double t) const;

    /// The Hull–White short-rate volatility on each ζ interval (size = number of nodes).
    std::vector<double> hullWhiteSigmas() const;

    /// P(t, T; x) / N(t, x) = D(T) exp(−H(T) x − ½ H(T)² ζ(t)), given D(T) and ζ(t).
    double deflatedBond(double dfT, double T, double zetaT, double x) const noexcept;
    /// N(t, x) given D(t) and ζ(t).
    double numeraire(double dfT, double t, double zetaT, double x) const noexcept;
    /// P(t, T; x) given D(t), D(T) and ζ(t).
    double bond(double dfT0, double t, double dfT1, double T, double zetaT, double x) const noexcept;

private:
    double kappa_;
    std::vector<double> times_;
    std::vector<double> zetas_;
};

} // namespace rl::models
