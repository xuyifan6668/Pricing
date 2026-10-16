#pragma once
// Per-node SABR calibration.
//
//   SabrNodeProblem                 the calibration problem for one (expiry, tenor): residuals, exact
//                                   Jacobian, weights and parameter bounds.
//   sabrInitialGuess, fitSabrSmile  starting point from the smile's skew and curvature; bounded
//                                   Levenberg–Marquardt; a verdict with failure flags.
//   smileQuotes, makeNodeProblem, fitCube   from the vol screen to a fitted SwaptionCube.
//
// Parameter vector z: with options.alphaFromAtm (the default) z = (ρ, ν) and α is re-solved from the ATM
// quote at every evaluation, so the fit reprices ATM exactly and has two unknowns; without it
// z = (α, ρ, ν). β and the shift are policy and never fitted: β and ρ both produce skew, and fitting
// both makes the problem ill-conditioned.
//
// Residuals are raw, in bp of normal vol: F_i(z) = 1e4 · (σ_SABR(K_i; z) − σ_i^mkt). Weights are kept
// separate: the report prints raw residuals, the solver minimises the weighted ones.

#include "rl/market/sabr_params.hpp"
#include "rl/market/swaption_cube.hpp"
#include "rl/market/vol_quotes.hpp"
#include "rl/math/calibration_problem.hpp"
#include "rl/math/lm.hpp"

#include <cstddef>
#include <string>
#include <utility>
#include <vector>

namespace rl::vol {

struct SmileQuote {
    double strike = 0.0;    ///< decimal
    double normalVol = 0.0; ///< decimal
};

enum class SmileWeighting {
    Equal,      ///< every quote counts the same in vol bp
    Vega,       ///< weight by Bachelier vega at the quote, relative to ATM (≈ premium error)
    InverseVol, ///< relative vol error
};

/// The bounds and the residual tolerance are part of the objective; their reasoning is above
/// fitSabrSmile in src/sabr_fit.cpp.
struct SabrFitOptions {
    double beta = 0.5;
    double shift = 0.0;
    bool alphaFromAtm = true;
    SmileWeighting weighting = SmileWeighting::Vega;
    double rhoMin = -0.9999;
    double rhoMax = 0.9999;
    double nuMin = 1e-4;
    double nuMax = 5.0;
    double alphaMin = 1e-8; ///< used only when !alphaFromAtm
    double alphaMax = 5.0;
    double maxResidualBp = 0.1; ///< a fit whose max |raw residual| exceeds this is not ok
    bool throwOnFailure = false; ///< CalibrationError (carrying the residuals in bp) instead of ok = false
    math::LmOptions lm;          ///< bounds are filled from the fields above by fitSabrSmile
};

struct SabrFitResult {
    market::SabrParams<double> params;
    std::vector<double> residualsBp; ///< raw, per quote, model − market
    double rmsBp = 0.0;
    double maxAbsBp = 0.0;
    math::LmStatus status = math::LmStatus::Failed;
    int iterations = 0;
    std::vector<std::string> flags; ///< e.g. "rho at bound", "max residual 0.31 bp > 0.1 bp"
    bool ok = false;
};

class SabrNodeProblem {
public:
    /// Throws InputError on: fewer quotes than parameters, duplicate strikes, non-positive or non-finite vols,
    /// expiry ≤ 0, (β > 0) a forward or strike with F + s ≤ 0 or K + s ≤ 0, or inconsistent options (bounds
    /// inverted or outside the parameter domain, β ∉ [0, 1], maxResidualBp ≤ 0). Quotes are stored sorted by strike.
    SabrNodeProblem(std::string id, std::vector<SmileQuote> quotes, double forward, double expiry, double atmNormalVol, SabrFitOptions options = {});

    const std::string& id() const noexcept { return id_; }
    const std::vector<SmileQuote>& quotes() const noexcept { return quotes_; }
    double forward() const noexcept { return forward_; }
    double expiry() const noexcept { return expiry_; }
    double atmNormalVol() const noexcept { return atm_; }
    const SabrFitOptions& options() const noexcept { return options_; }

    // --- math::CalibrationProblem ---------------------------------------------------------------------
    std::size_t parameterCount() const noexcept { return options_.alphaFromAtm ? 2 : 3; }
    std::size_t residualCount() const noexcept { return quotes_.size(); }
    /// "ATM", "-100bp", "+25bp", … by offset from the forward.
    const std::vector<std::string>& instrumentIds() const noexcept { return ids_; }
    /// Raw residuals in bp of normal vol. Throws NumericalError if α cannot be re-solved from ATM.
    math::Vector residuals(const math::Vector& z) const;
    /// Per-residual weights S_i (— the weighting; see src/sabr_fit.cpp).
    math::Vector scaling(const math::Vector& z) const;

    // --- helpers ---------------------------------------------------------------------------------------
    /// Full parameters from z (α re-solved when alphaFromAtm). Throws InputError for a z of the wrong size.
    market::SabrParams<double> params(const math::Vector& z) const;
    /// z from full parameters (drops α when alphaFromAtm).
    math::Vector toVector(const market::SabrParams<double>& p) const;
    /// ∂residuals/∂z, residualCount × parameterCount, exact (Dual through the formula and the ATM constraint).
    math::Matrix jacobian(const math::Vector& z) const;
    /// (lower, upper) bounds on z from the options, in z's layout.
    std::pair<math::Vector, math::Vector> bounds() const;

private:
    std::string id_;
    std::vector<SmileQuote> quotes_;
    double forward_;
    double expiry_;
    double atm_;
    SabrFitOptions options_;
    std::vector<std::string> ids_;
};

/// a starting z inside bounds() from the quotes alone.
math::Vector sabrInitialGuess(const SabrNodeProblem& problem);

/// calibrate one node.
SabrFitResult fitSabrSmile(const SabrNodeProblem& problem);

// --- plumbing -----------------------------------------------------------------------------------

/// The node's ATM quote (strike = forward) and smile quotes (strike = forward + offset), sorted by strike.
std::vector<SmileQuote> smileQuotes(const market::SwaptionVolNodeQuote& node);

SabrNodeProblem makeNodeProblem(const market::SwaptionVolNodeQuote& node, const SabrFitOptions& options);

struct CubeFit {
    market::SwaptionCube<double> cube;
    std::vector<std::string> nodeIds;  ///< row-major (expiry, tenor), aligned with `fits`
    std::vector<SabrFitResult> fits;
};

/// Fit every node of a complete vol grid and assemble the cube (row-major by the quotes' sorted axes). Nodes
/// that are not ok still enter the cube with their best parameters; the caller reads `fits` for flags.
CubeFit fitCube(const market::SwaptionVolQuotes& quotes, const SabrFitOptions& options);

} // namespace rl::vol
