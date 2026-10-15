#pragma once
// Shifted SABR: Hagan's normal-vol expansion, its lognormal twin, the ATM constraint and the gradient.
//
// Model (market/sabr_params.hpp): d(F+s) = α_t (F+s)^β dW, dα_t = ν α_t dZ, dW dZ = ρ dt.
//
// Normal-vol form (primary) — Hagan, Kumar, Lesniewski & Woodward (2002), "Managing Smile Risk", Wilmott
// Magazine, Appendix B, eq. (B.69):
//
//   σ_N(K) = α (1−β)(f − k) / (f^{1−β} − k^{1−β}) · ζ / x̂(ζ)
//            · { 1 + [ −β(2−β) α² / (24 f_av^{2−2β}) + ρ α ν β / (4 f_av^{1−β}) + (2 − 3ρ²) ν² / 24 ] T }
//
//   f = F + s,  k = K + s,  f_av = √(f k),
//   ζ = (ν/α) (f^{1−β} − k^{1−β}) / (1−β)            (→ (ν/α) ln(f/k) at β = 1; (ν/α)(F − K) at β = 0)
//   x̂(ζ) = ln( (√(1 − 2ρζ + ζ²) + ζ − ρ) / (1 − ρ) )
//
// DESIGN: ζ is written in its integral form ζ = (ν/α)∫_k^f u^{−β} du. The paper's (B.69b) writes the
// leading-order equivalent (ν/α)(f − k)/f_av^β; the two agree to the order of the expansion and both are in
// production use. The integral form is exact at β = 0 and β = 1, and it is what the ν → 0 limit of the first
// factor is built from, so the formula is internally consistent. Say this at the whiteboard.
//
// Lognormal form — ibid. eqs. (2.17a,b), for the notebook's cross-check only:
//   σ_B(K) = α / { (fk)^{(1−β)/2} [1 + (1−β)² ln²(f/k)/24 + (1−β)⁴ ln⁴(f/k)/1920] } · ζ/x̂(ζ)
//            · { 1 + [ (1−β)² α² / (24 (fk)^{1−β}) + ρβνα / (4 (fk)^{(1−β)/2}) + (2−3ρ²)ν²/24 ] T },
//   ζ = (ν/α)(fk)^{(1−β)/2} ln(f/k).
//
// Limits handled without a branch that changes the derivative: ζ → 0 (series ζ/x̂ = 1 − ρζ/2 + (2−3ρ²)ζ²/12),
// f → k (series of the first factor), β = 0 (no logs or powers, so F + s and K + s may be ≤ 0) and β = 1.
// Branches are on value() only, and every series branch has the right first derivative.
//
// Where it fails (the honest part, notebook 07 and vol/arbitrage.hpp): the expansion is in ε² ~ α²T f^{2β−2}
// and ν²T. For long expiries with large ν, and at low strikes when β > 0 (f_av^{2−2β} → 0 drives the α² term
// to −∞), σ_N(K) can make the implied density negative and even the vol itself negative. The production fix is
// Hagan, Kumar, Lesniewski & Woodward (2014), "Arbitrage-free SABR", Wilmott Magazine; v1 reports, does not repair.
//
// Active set: α, ρ, ν and the forward are T; strike, expiry, β and shift are passive.

#include "rl/market/sabr_params.hpp"
#include "rl/math/scalar.hpp"

namespace rl::vol {

/// Hagan normal vol (B.69) at strike K. Throws InputError on invalid parameters, a negative expiry, or (β > 0)
/// F + s ≤ 0 or K + s ≤ 0. May return a non-positive number deep in the wing for β > 0 — that is the
/// expansion failing, reported by vol/arbitrage.hpp, not an error here.
template <math::Scalar T>
T sabrNormalVol(const market::SabrParams<T>& params, const T& forward, double strike, double expiry);

/// Hagan lognormal (shifted Black) vol (2.17a). Requires F + s > 0 and K + s > 0 for every β.
template <math::Scalar T>
T sabrLognormalVol(const market::SabrParams<T>& params, const T& forward, double strike, double expiry);

/// σ_N at K = F in closed form: α f^β {1 + [−β(2−β)α²/(24 f^{2−2β}) + ραβν/(4 f^{1−β}) + (2−3ρ²)ν²/24] T},
/// a cubic in α. Equals sabrNormalVol(params, F, value(F), T) to rounding.
template <math::Scalar T>
T sabrAtmNormalVol(const market::SabrParams<T>& params, const T& forward, double expiry);

/// α such that sabrAtmNormalVol = atmNormalVol: the smallest positive root of the ATM cubic (West 2005,
/// "Calibration of the SABR model in illiquid markets", Applied Mathematical Finance 12(4), §2), which is
/// bracketed by 0 and the cubic's first positive stationary point.
///
/// DESIGN: the root is found on double (Newton with bisection safeguard), then ONE Newton step is taken
/// in T arithmetic at the converged root, α_T = α₀ − g(α₀; θ_T)/g'(α₀; θ_T). The value moves by ~1e-17; the
/// derivative recorded is exactly the implicit-function-theorem derivative dα/dθ = −(∂g/∂θ)/(∂g/∂α), so vega
/// through the ATM constraint is right under aad::Real today and under Dual for the fit's Jacobian, with no
/// taping of the solver. the Jamshidian root uses the same pattern.
///
/// Throws InputError on invalid inputs; NumericalError if the ATM vol exceeds the maximum of the cubic (the
/// expansion cannot produce it for these ρ, ν, β — a stale or mistyped ATM quote, or a very long expiry).
template <math::Scalar T>
T alphaFromAtmNormalVol(const T& atmNormalVol, const T& forward, double expiry, double beta, const T& rho, const T& nu, double shift);

template <math::Scalar T>
struct SabrVolGradient {
    T vol;
    T dAlpha;
    T dRho;
    T dNu;
};

/// σ_N(K) and ∂σ_N/∂(α, ρ, ν) at fixed forward. Exact: forward-mode Dual<T, 3> through sabrNormalVol itself
/// (math/dual.hpp), so the gradient cannot drift from the formula.
template <math::Scalar T>
SabrVolGradient<T> sabrNormalVolGradient(const market::SabrParams<T>& params, const T& forward, double strike, double expiry);

template <math::Scalar T>
struct SabrAtmConstrainedGradient {
    T vol;
    T dRho;   ///< total derivative, α re-solved from the ATM vol
    T dNu;
    T alpha;
    T dAlphaDRho;
    T dAlphaDNu;
};

/// σ_N(K) and its total derivatives in (ρ, ν) when α is tied to a fixed ATM normal vol by alphaFromAtmNormalVol.
/// The Jacobian a two-parameter SABR fit needs (vol/sabr_fit.hpp). Forward-mode Dual<T, 2>.
template <math::Scalar T>
SabrAtmConstrainedGradient<T> sabrNormalVolGradientAtmConstrained(const T& atmNormalVol, const T& forward, double strike, double expiry, double beta,
                                                                  const T& rho, const T& nu, double shift);

} // namespace rl::vol
