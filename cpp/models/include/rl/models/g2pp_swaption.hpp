#pragma once
// European options on a sum of bonds under G2++.
//
// The deflated payoff of exercising at t into flows {a_k at T_k} is
//
//     g(x_1, x_2) = Σ_k c_k exp(−L_k · x),     c_k = a_k D(T_k) exp(−½ L_kᵀ C(t) L_k),
//
// with x ~ N(0, C(t)), and the option is worth E[g⁺]. The loadings are L_k = H(T_k) − H(t): the
// payoff is written under the numeraire P(·, t), where every term is well scaled (g2pp.hpp). With two states there is no single exercise
// point, so Jamshidian's decomposition does not apply directly. It does apply one dimension down:
// pick one direction y in the state space, condition on what is independent of it, and for a swap g
// crosses zero once in y. The inner expectation is then a closed-form sum of normal cdfs, and what
// is left is a one-dimensional integral of a smooth function, done by Gauss–Hermite quadrature.
// Brigo & Mercurio (eq. 4.31) take y = x_2; the direction used here is chosen so the integral stays
// accurate when the factors are highly correlated, and is exact in the one-factor limit
// (σ_2 = 0 or |ρ| = 1), where the price is the one-dimensional closed form.

#include "rl/models/flows.hpp"
#include "rl/models/g2pp.hpp"

#include <cstdint>
#include <vector>

namespace rl::models {

/// E[g(x)⁺] at exercise time `t`: the price today of the right to receive `flows` then. The payoff
/// must be swap-like (one sign change in x_2 for each x_1). `quadraturePoints` is the Gauss–Hermite
/// order of the outer integral.
double g2ppOptionValue(const G2pp& model, double t, const std::vector<BondFlow>& flows, int quadraturePoints = 24);

/// The same expectation by Simpson's rule on a two-dimensional grid over ±8 standard deviations.
/// Makes no assumption about the payoff's shape; used to check the semi-closed form.
double g2ppOptionValueByIntegration(const G2pp& model, double t, const std::vector<BondFlow>& flows, int pointsPerAxis = 801);

/// Deflated exercise value g(x_1, x_2) at time t.
double g2ppDeflatedValue(const G2pp& model, const G2pp::Covariance& c, const std::vector<BondFlow>& flows, double x1, double x2) noexcept;

/// Model price of the European swaption with this exercise right.
double g2ppSwaptionPrice(const G2pp& model, const ExerciseRight& right);

/// Price of the option to receive 1 at `expiry` and pay `strikeAmount` at `maturity` from the
/// two-factor bond-option formula: Black's formula on the forward bond with variance
/// (H(maturity) − H(expiry))ᵀ C(expiry) (H(maturity) − H(expiry)). Exact and closed form.
double g2ppZeroBondPut(const G2pp& model, double expiry, double maturity, double dfExpiry, double dfMaturity, double strikeAmount);

} // namespace rl::models
