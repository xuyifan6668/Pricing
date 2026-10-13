#pragma once
// Shifted Black-76 (shifted lognormal) model for an option on a forward rate.
//
//   d(F + s) = σ (F + s) dW;   f = F + s,  k = K + s,  d₁ = [ln(f/k) + σ²T/2] / (σ√T),  d₂ = d₁ − σ√T
//   C = f Φ(d₁) − k Φ(d₂),     P = k Φ(−d₂) − f Φ(−d₁)
//   Δ_C = Φ(d₁),  Δ_P = −Φ(−d₁),  vega = f φ(d₁) √T,  Γ = φ(d₁) / (f σ √T)
//
// Source: Black (1976), "The pricing of commodity contracts"; the shift is the displaced-diffusion model of
// Rubinstein (1983). Reference value (tested): Haug, "The Complete Guide to Option Pricing Formulas" 2e (2007),
// §1.1.4: F = 19, K = 19, T = 0.75, r = 10%, σ = 28% gives a discounted call and put of 1.7011.
//
// Undiscounted, as bachelier.hpp. The shift s is passive policy. f = F + s must be positive (the model
// does not exist otherwise: InputError). k = K + s ≤ 0 is legitimate and vol-independent: the call is worth
// F − K and the put nothing, because f > 0 ≥ k almost surely.

#include "rl/math/scalar.hpp"
#include "rl/vol/bachelier.hpp" // ImpliedVolOptions
#include "rl/vol/option_type.hpp"

namespace rl::vol {

template <math::Scalar T>
T blackPrice(const T& forward, double strike, double expiry, const T& lognormalVol, double shift, OptionType type);

template <math::Scalar T>
T blackDelta(const T& forward, double strike, double expiry, const T& lognormalVol, double shift, OptionType type);

/// ∂price/∂σ, identical for calls and puts; zero when K + s ≤ 0.
template <math::Scalar T>
T blackVega(const T& forward, double strike, double expiry, const T& lognormalVol, double shift);

template <math::Scalar T>
T blackGamma(const T& forward, double strike, double expiry, const T& lognormalVol, double shift);

/// σ such that blackPrice(F, K, T, σ, s, type) = price. Solves on the out-of-the-money option (by parity the
/// time value is the same) with safeguarded Newton on log(price) in σ inside an expanded bracket, falling back
/// to Brent on the price itself. Returns 0 at intrinsic. Throws InputError if expiry ≤ 0, F + s ≤ 0, K + s ≤ 0
/// (the price carries no vol information), the price is below intrinsic, or at/above the no-arbitrage bound
/// (min(F + s, K + s) for the OTM option).
double blackImpliedVol(double price, double forward, double strike, double expiry, double shift, OptionType type, const ImpliedVolOptions& options = {});

} // namespace rl::vol
