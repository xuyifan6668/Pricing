#pragma once
// Static-arbitrage diagnostic for a SABR smile (double only). Reports; never repairs (v1).
//
// Breeden & Litzenberger (1978): with C(K) the undiscounted call on the forward, the risk-neutral density is
// p(K) = ∂²C/∂K², and absence of static arbitrage in strike requires
//   p(K) ≥ 0            (butterflies have non-negative value)
//   −1 ≤ ∂C/∂K ≤ 0      (call spreads are worth between 0 and the strike width)
// Here C(K) = bachelierPrice(F, K, T, σ_SABR(K)) and both derivatives are central differences with step h.
//
// DESIGN: a diagnostic, not a fix. The production remedy is the arbitrage-free SABR of Hagan, Kumar, Lesniewski &
// Woodward (2014), "Arbitrage-free SABR", Wilmott Magazine (solve the effective 1-D forward PDE for the density),
// or an extrapolation that replaces the wing beyond the first bad strike. Either changes every wing price; v1's
// job is to say where the Hagan expansion stops being usable, so that nobody replicates a CMS off it unknowingly.
//
// What gets reported and why the scan starts at the forward: for β > 0 and no shift, Hagan's density goes negative
// as K → 0 on every long-expiry node (the α² correction term diverges as (FK)^{β−1}), so "is there a negative
// density somewhere" is always yes. The number a desk uses is the highest strike below the forward where the
// density first goes negative — the edge of the usable smile — and it moves up sharply with ν²T.
//
// Numerical choices: h = 0.5 bp by default. The second difference carries rounding noise ≈ 4·eps·F/h² ≈ 1e-8 and
// truncation ≈ h² p''/12 (≈ 1e-5 relative on a 1M smile); the densities that matter are O(1)–O(100) per unit
// rate, so a negativity tolerance of 1e-6 separates noise from arbitrage by several orders of magnitude.

#include "rl/market/sabr_params.hpp"

#include <cstddef>
#include <optional>
#include <utility>
#include <vector>

namespace rl::vol {

struct DensityPoint {
    double strike = 0.0;
    double normalVol = 0.0; ///< σ_SABR(K); may be ≤ 0 where the expansion has broken down
    double callPrice = 0.0; ///< NaN when the vol at K − h, K or K + h is ≤ 0
    double dCalldK = 0.0;   ///< NaN as above
    double density = 0.0;   ///< NaN as above
};

/// C, ∂C/∂K and ∂²C/∂K² on the given strikes. Throws InputError if expiry ≤ 0, h ≤ 0, or (β > 0) a strike with
/// K − h + shift ≤ 0.
std::vector<DensityPoint> riskNeutralDensity(const market::SabrParams<double>& params, double forward, double expiry, const std::vector<double>& strikes,
                                             double h = 5e-5);

struct ArbitrageCheck {
    double kMin = 0.0;
    double kMax = 0.0;
    std::size_t points = 0;
    double minDensity = 0.0;
    double minDensityStrike = 0.0;
    /// Scanning DOWN from the forward: the first strike with density < −tol or an unusable vol (≤ 0).
    std::optional<double> firstNegativeDensityBelowForward;
    /// Scanning UP from the forward: the first strike with density < −tol or an unusable vol.
    std::optional<double> firstNegativeDensityAboveForward;
    /// Any strike where ∂C/∂K leaves [−1, 0] by more than 1e-10.
    bool callSpreadViolated = false;
    /// Any strike where σ_SABR ≤ 0.
    bool nonPositiveVol = false;
    /// Trapezoid ∫p dK and ∫K p dK over the grid (≈ 1 and ≈ F on a grid that covers the distribution).
    double integratedMass = 0.0;
    double impliedMean = 0.0;

    bool arbitrageFree() const noexcept {
        return !firstNegativeDensityBelowForward && !firstNegativeDensityAboveForward && !callSpreadViolated && !nonPositiveVol;
    }
};

/// Summarise a density grid (strikes ascending) around `forward`.
ArbitrageCheck checkDensity(const std::vector<DensityPoint>& density, double forward, double densityTol = 1e-6);

/// riskNeutralDensity on a uniform grid of n strikes over [kMin, kMax], then checkDensity.
ArbitrageCheck checkSabrSmile(const market::SabrParams<double>& params, double forward, double expiry, double kMin, double kMax, std::size_t n);

/// n uniformly spaced strikes from kMin to kMax inclusive (n ≥ 2, kMax > kMin).
std::vector<double> strikeGrid(double kMin, double kMax, std::size_t n);

/// The default diagnostic range F ± nStdDev · σ_ATM √T, with the lower end floored at −shift + floor when β > 0
/// (the formula is undefined at K + s ≤ 0).
std::pair<double, double> diagnosticStrikeRange(const market::SabrParams<double>& params, double forward, double expiry, double atmNormalVol,
                                                double nStdDev = 4.0, double floor = 1e-4);

} // namespace rl::vol
