#pragma once
// Bachelier (normal) model for an option on a forward rate.
//
//   dF = σ_N dW  under the annuity measure;   d = (F − K) / (σ_N √T)
//   C = (F − K) Φ(d) + σ_N √T φ(d),     P = (K − F) Φ(−d) + σ_N √T φ(d)
//   Δ_C = Φ(d),  Δ_P = −Φ(−d),  ν = √T φ(d),  Γ = φ(d) / (σ_N √T)
//
// Source: Bachelier (1900); in rates form Andersen & Piterbarg, "Interest Rate Modeling" vol. 1 (2010), §7.1.1.
// Undiscounted: multiply by the annuity for a swaption premium. F and K may be negative (that is why rates
// desks quote normal vols). ATM check: C = P = σ_N √T / √(2π).
//
// Active set: the forward and the vol are T (the forward comes off the curve, the vol off the cube); strike and
// expiry are trade terms and passive. The implied-vol solver is double-only: nothing under T = Real
// inverts a price, because the cube stores normal vols and SABR returns one directly.
//
// Numerical note: the OTM time value s φ(d) − |F − K| Φ(−|d|) cancels relatively as |d| grows (a factor ≈ d²),
// but the vol's sensitivity to it grows by the same factor, so the implied vol keeps ≈ machine precision until
// φ(d) underflows near |d| ≈ 38. The tests round-trip to 1e-12 at ±6 standard deviations.

#include "rl/math/scalar.hpp"
#include "rl/vol/option_type.hpp"

namespace rl::vol {

template <math::Scalar T>
T bachelierPrice(const T& forward, double strike, double expiry, const T& normalVol, OptionType type);

/// ∂price/∂F.
template <math::Scalar T>
T bachelierDelta(const T& forward, double strike, double expiry, const T& normalVol, OptionType type);

/// ∂price/∂σ_N (identical for calls and puts).
template <math::Scalar T>
T bachelierVega(const T& forward, double strike, double expiry, const T& normalVol);

/// ∂²price/∂F² (identical for calls and puts). Zero when σ√T = 0 and F ≠ K (the Dirac at F = K is not represented).
template <math::Scalar T>
T bachelierGamma(const T& forward, double strike, double expiry, const T& normalVol);

struct ImpliedVolOptions {
    double relTol = 1e-14; ///< stop when the step in σ√T is below relTol × the bracket scale
    int maxIter = 200;
};

/// σ_N such that bachelierPrice(F, K, T, σ_N, type) = price. Solves on the out-of-the-money time value
/// (price − intrinsic, equal for the call and the put by parity) with safeguarded Newton on log(time value)
/// inside the exact bracket  s/√(2π) − |F−K|/2 ≤ h(s) ≤ s/√(2π),  s = σ√T.
/// Returns 0 when the price equals intrinsic. Throws InputError if expiry ≤ 0, inputs are not finite, or the
/// price is below intrinsic; NumericalError if the solver fails (it should not).
double bachelierImpliedVol(double price, double forward, double strike, double expiry, OptionType type, const ImpliedVolOptions& options = {});

} // namespace rl::vol
