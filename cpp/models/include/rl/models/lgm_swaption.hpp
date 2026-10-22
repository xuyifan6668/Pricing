#pragma once
// European options on a sum of bonds under LGM, and the same swaption on the market's terms.
//
// Under the model, the deflated payoff of exercising at t into flows {a_k at T_k} is
//
//     g(x) = Σ_k a_k D(T_k) exp(−H(T_k) x − ½ H(T_k)² ζ(t)),        x ~ N(0, ζ(t)),
//
// and the option is worth E[g(x)⁺]. For a swap, g crosses zero exactly once (the swap is worth
// more to the payer the higher rates are), at x*. Then the option is a sum of bond options that are
// all exercised together, each with a closed form (Jamshidian's decomposition):
//
//     V = Σ_k a_k D(T_k) Φ( ∓ (x* + H(T_k) ζ) / √ζ ),      − if g > 0 above x*, + if below.
//
// On the market's terms the same swaption is quoted as a normal vol σ_N of the forward swap rate:
//     V = A(0) · Bachelier(F, K, t, σ_N),     A = annuity, F = forward par rate.
// The two meet in calibration, where ζ(t) is chosen so the model price equals the market price.

#include "rl/models/flows.hpp"
#include "rl/models/lgm.hpp"

#include <vector>

namespace rl::models {

/// E[(Σ a_k D_k exp(−H_k x − ½ H_k² ζ))⁺] for x ~ N(0, ζ): the price today of the right to receive
/// `flows` at the time whose state variance is `zeta`. Uses the single-crossing closed form; if the
/// payoff does not change sign (always or never worth exercising) returns max(Σ a_k D_k, 0).
/// Throws NumericalError if the payoff crosses zero more than once (not a swap-like payoff).
double lgmOptionValue(const Lgm1F& model, double zeta, const std::vector<BondFlow>& flows);

/// The same expectation by direct numerical integration (composite Simpson over ±10 standard
/// deviations). Makes no assumption about the payoff's shape; used to check the closed form.
double lgmOptionValueByIntegration(const Lgm1F& model, double zeta, const std::vector<BondFlow>& flows, int points = 4001);

/// Deflated exercise value g(x) at state x.
double deflatedValue(const Lgm1F& model, double zeta, const std::vector<BondFlow>& flows, double x) noexcept;

/// The state x* at which exercising is worth exactly zero, and which side is in the money.
struct ExerciseBoundary {
    bool exists = false;    ///< false: the payoff has one sign for every plausible state
    double x = 0.0;
    bool exerciseAbove = true;
};
ExerciseBoundary exerciseBoundary(const Lgm1F& model, double zeta, const std::vector<BondFlow>& flows);

/// Model price of the European swaption with this exercise right.
double lgmSwaptionPrice(const Lgm1F& model, const ExerciseRight& right);

/// Market price from a normal vol: annuity × Bachelier on the forward par rate.
double marketSwaptionPrice(const ExerciseRight& right, double normalVol);

/// The normal vol at which the market formula gives `price` (0 if the price is intrinsic).
double impliedNormalVol(const ExerciseRight& right, double price);

} // namespace rl::models
