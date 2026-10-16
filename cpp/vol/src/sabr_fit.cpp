#include "rl/vol/sabr_fit.hpp"

#include "rl/core/error.hpp"
#include "rl/vol/bachelier.hpp"
#include "rl/vol/sabr.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>

namespace rl::vol {

// --- SabrNodeProblem: construction and residuals/Jacobian ---------------------------------------

namespace {
void validateOptions(const SabrFitOptions& o) {
    RL_REQUIRE(std::isfinite(o.beta) && o.beta >= 0.0 && o.beta <= 1.0, "SabrFitOptions: beta must lie in [0, 1]");
    RL_REQUIRE(std::isfinite(o.shift), "SabrFitOptions: shift must be finite");
    RL_REQUIRE(o.rhoMin > -1.0 && o.rhoMin < o.rhoMax && o.rhoMax < 1.0, "SabrFitOptions: need -1 < rhoMin < rhoMax < 1");
    RL_REQUIRE(o.nuMin >= 0.0 && o.nuMin < o.nuMax && std::isfinite(o.nuMax), "SabrFitOptions: need 0 <= nuMin < nuMax");
    RL_REQUIRE(o.alphaMin > 0.0 && o.alphaMin < o.alphaMax && std::isfinite(o.alphaMax), "SabrFitOptions: need 0 < alphaMin < alphaMax");
    RL_REQUIRE(o.maxResidualBp > 0.0, "SabrFitOptions: maxResidualBp must be positive");
}

std::string offsetLabel(double strike, double forward) {
    const double bp = (strike - forward) * 1e4;
    if (std::fabs(bp) < 0.5) return "ATM";
    const long long r = std::llround(bp);
    return (r > 0 ? "+" : "") + std::to_string(r) + "bp";
}
} // namespace

SabrNodeProblem::SabrNodeProblem(std::string id, std::vector<SmileQuote> quotes, double forward, double expiry, double atmNormalVol, SabrFitOptions options)
    : id_(std::move(id)), quotes_(std::move(quotes)), forward_(forward), expiry_(expiry), atm_(atmNormalVol), options_(std::move(options)) {
    withContext("SABR node " + id_, [&] {
        validateOptions(options_);
        RL_REQUIRE(std::isfinite(forward_), "forward must be finite");
        RL_REQUIRE(std::isfinite(expiry_) && expiry_ > 0.0, "expiry must be positive");
        RL_REQUIRE(std::isfinite(atm_) && atm_ > 0.0, "ATM normal vol must be positive");
        if (options_.beta > 0.0) RL_REQUIRE(forward_ + options_.shift > 0.0, "forward + shift must be positive when beta > 0");
        std::sort(quotes_.begin(), quotes_.end(), [](const SmileQuote& a, const SmileQuote& b) { return a.strike < b.strike; });
        for (std::size_t i = 0; i < quotes_.size(); ++i) {
            const auto& q = quotes_[i];
            RL_REQUIRE(std::isfinite(q.strike) && std::isfinite(q.normalVol) && q.normalVol > 0.0, "quotes need finite strikes and positive vols");
            if (options_.beta > 0.0) RL_REQUIRE(q.strike + options_.shift > 0.0, "strike + shift must be positive when beta > 0");
            if (i > 0) RL_REQUIRE(q.strike > quotes_[i - 1].strike, "duplicate strike");
        }
        RL_REQUIRE(quotes_.size() >= parameterCount(), "need at least as many quotes as parameters (" + std::to_string(parameterCount()) + ")");
    });
    ids_.reserve(quotes_.size());
    for (const auto& q : quotes_) ids_.push_back(offsetLabel(q.strike, forward_));
}

market::SabrParams<double> SabrNodeProblem::params(const math::Vector& z) const {
    RL_REQUIRE(static_cast<std::size_t>(z.size()) == parameterCount(), "SABR node " + id_ + ": parameter vector has the wrong size");
    market::SabrParams<double> p;
    p.beta = options_.beta;
    p.shift = options_.shift;
    if (options_.alphaFromAtm) {
        p.rho = z[0];
        p.nu = z[1];
        p.alpha = alphaFromAtmNormalVol<double>(atm_, forward_, expiry_, p.beta, p.rho, p.nu, p.shift);
    } else {
        p.alpha = z[0];
        p.rho = z[1];
        p.nu = z[2];
    }
    return p;
}

math::Vector SabrNodeProblem::toVector(const market::SabrParams<double>& p) const {
    math::Vector z(static_cast<Eigen::Index>(parameterCount()));
    if (options_.alphaFromAtm) z << p.rho, p.nu;
    else z << p.alpha, p.rho, p.nu;
    return z;
}

math::Vector SabrNodeProblem::residuals(const math::Vector& z) const {
    const market::SabrParams<double> p = params(z);
    math::Vector r(static_cast<Eigen::Index>(quotes_.size()));
    for (std::size_t i = 0; i < quotes_.size(); ++i)
        r[static_cast<Eigen::Index>(i)] = 1e4 * (sabrNormalVol<double>(p, forward_, quotes_[i].strike, expiry_) - quotes_[i].normalVol);
    return r;
}

math::Matrix SabrNodeProblem::jacobian(const math::Vector& z) const {
    const market::SabrParams<double> p = params(z);
    const Eigen::Index m = static_cast<Eigen::Index>(quotes_.size());
    math::Matrix J(m, static_cast<Eigen::Index>(parameterCount()));
    for (Eigen::Index i = 0; i < m; ++i) {
        const double k = quotes_[static_cast<std::size_t>(i)].strike;
        if (options_.alphaFromAtm) {
            const auto g = sabrNormalVolGradientAtmConstrained<double>(atm_, forward_, k, expiry_, p.beta, p.rho, p.nu, p.shift);
            J(i, 0) = 1e4 * g.dRho;
            J(i, 1) = 1e4 * g.dNu;
        } else {
            const auto g = sabrNormalVolGradient<double>(p, forward_, k, expiry_);
            J(i, 0) = 1e4 * g.dAlpha;
            J(i, 1) = 1e4 * g.dRho;
            J(i, 2) = 1e4 * g.dNu;
        }
    }
    return J;
}

std::pair<math::Vector, math::Vector> SabrNodeProblem::bounds() const {
    const auto n = static_cast<Eigen::Index>(parameterCount());
    math::Vector lo(n), hi(n);
    if (options_.alphaFromAtm) {
        lo << options_.rhoMin, options_.nuMin;
        hi << options_.rhoMax, options_.nuMax;
    } else {
        lo << options_.alphaMin, options_.rhoMin, options_.nuMin;
        hi << options_.alphaMax, options_.rhoMax, options_.nuMax;
    }
    return {lo, hi};
}

// --- the calibration objective ------------------------------------------------------------------------------

namespace {
/// Vega weights are floored here. At ±100 bp on a 1M expiry (σ√T ≈ 30 bp, so |d| ≈ 3.3) the raw vega
/// weight is exp(−d²/2) ≈ 0.004: the fit would ignore the wings entirely and ν would be set by the
/// ATM curvature alone. A 10% floor keeps every quoted strike in the objective.
constexpr double kVegaWeightFloor = 0.10;
/// A parameter this close to a bound (relative to the bound's width) is "on the bound".
constexpr double kBoundTolerance = 1e-6;
/// Residual returned when α cannot be re-solved from the ATM vol at a trial point, so LM rejects the step.
constexpr double kPenaltyBp = 1e3;
} // namespace

// Weights S; the fit minimises ½‖S ∘ F‖² with F the raw residuals in bp of normal vol.
//   Equal       1 everywhere: every quote counts the same in vol.
//   Vega        Bachelier vega at the quote relative to ATM, floored: a vol error of x bp at strike K
//               costs vega(K)·x of premium, so this weights (roughly) premium error. Right for a book
//               whose risk sits near the money.
//   InverseVol  σ_ATM/σ_i: relative vol error.
// DESIGN: the weights are evaluated at the MARKET vols and do not depend on z. Weights that moved with
// the model would change the objective under the solver (its Gauss–Newton model ignores ∂S/∂z) and can
// create spurious minima.
math::Vector SabrNodeProblem::scaling(const math::Vector&) const {
    const auto m = static_cast<Eigen::Index>(quotes_.size());
    math::Vector s = math::Vector::Ones(m);
    if (options_.weighting == SmileWeighting::Equal) return s;
    const double atmVega = bachelierVega<double>(forward_, forward_, expiry_, atm_);
    for (Eigen::Index i = 0; i < m; ++i) {
        const SmileQuote& q = quotes_[static_cast<std::size_t>(i)];
        if (options_.weighting == SmileWeighting::Vega)
            s[i] = std::max(kVegaWeightFloor, bachelierVega<double>(forward_, q.strike, expiry_, q.normalVol) / atmVega);
        else
            s[i] = atm_ / q.normalVol;
    }
    return s;
}

// Starting point from the shape of the smile at the money.
//
// Expanding Hagan's normal-vol formula in x = K − F (with σ₀ the ATM vol and F̃ = F + shift):
//
//   σ_N(K) ≈ σ₀ [ 1 + ( β/(2F̃) + ρν/(2σ₀) ) x + ... ]
//
//   skew       s = ∂σ/∂K   ≈ σ₀β/(2F̃)  +  ρν/2
//   curvature  c = ∂²σ/∂K² ≈ (2 − 3ρ²)ν²/(6σ₀)  +  βρν/(2F̃)  +  σ₀β(β−1)/(4F̃²)
//
// The first term of the skew is the BACKBONE: with β > 0 the local vol is higher at higher rates even
// with no correlation. What is left after removing it is ρν/2: correlation between the rate and its
// vol tilts the smile. The curvature is driven by ν² (vol of vol fattens both tails). So
//
//   A := ρν       = 2 (s − σ₀β/(2F̃))
//   B := (2−3ρ²)ν² = 6σ₀ (c − βA/(2F̃) − σ₀β(β−1)/(4F̃²))        ⇒   ν² = (B + 3A²)/2,   ρ = A/ν.
//
// s and c come from a quadratic least-squares fit through all the quotes, which is less noisy than
// differencing the two strikes next to the money.
math::Vector sabrInitialGuess(const SabrNodeProblem& problem) {
    const auto& q = problem.quotes();
    const SabrFitOptions& o = problem.options();
    const double sigma0 = problem.atmNormalVol();
    const double fShift = problem.forward() + o.shift;

    double skew = 0.0, curvature = 0.0;
    bool haveCurvature = q.size() >= 3;
    if (haveCurvature) {
        // σ(x) ≈ a0 + a1 x + a2 x², x in units of 1% so the normal equations are well scaled.
        Eigen::MatrixXd X(static_cast<Eigen::Index>(q.size()), 3);
        Eigen::VectorXd y(static_cast<Eigen::Index>(q.size()));
        for (std::size_t i = 0; i < q.size(); ++i) {
            const double x = (q[i].strike - problem.forward()) * 100.0;
            X.row(static_cast<Eigen::Index>(i)) << 1.0, x, x * x;
            y[static_cast<Eigen::Index>(i)] = q[i].normalVol;
        }
        const Eigen::Vector3d a = X.colPivHouseholderQr().solve(y);
        skew = a[1] * 100.0;
        curvature = 2.0 * a[2] * 1e4;
    } else {
        skew = (q.back().normalVol - q.front().normalVol) / (q.back().strike - q.front().strike);
    }

    const double backboneSkew = o.beta > 0.0 ? sigma0 * o.beta / (2.0 * fShift) : 0.0;
    const double A = 2.0 * (skew - backboneSkew);
    double nu = 0.3; // two quotes carry no curvature information: a typical rates vol-of-vol
    if (haveCurvature) {
        const double backboneCurv = o.beta > 0.0 ? o.beta * A / (2.0 * fShift) + sigma0 * o.beta * (o.beta - 1.0) / (4.0 * fShift * fShift) : 0.0;
        const double B = 6.0 * sigma0 * (curvature - backboneCurv);
        const double nu2 = 0.5 * (B + 3.0 * A * A);
        // A concave smile gives ν² ≤ 0: SABR cannot produce it. Start low and let the fit flag the node.
        nu = nu2 > 0.0 ? std::sqrt(nu2) : 0.05;
    }
    const auto [lo, hi] = problem.bounds();
    const Eigen::Index iRho = o.alphaFromAtm ? 0 : 1, iNu = iRho + 1;
    // Stay strictly inside the box: a start on a bound stalls a clamping solver.
    nu = std::clamp(nu, lo[iNu] + 0.01 * (hi[iNu] - lo[iNu]), lo[iNu] + 0.5 * (hi[iNu] - lo[iNu]));
    const double rho = std::clamp(A / nu, std::max(-0.9, lo[iRho] + 1e-3), std::min(0.9, hi[iRho] - 1e-3));

    math::Vector z(static_cast<Eigen::Index>(problem.parameterCount()));
    if (o.alphaFromAtm) {
        z << rho, nu;
    } else {
        double alpha = sigma0 / std::pow(fShift > 0.0 ? fShift : 1.0, o.beta); // backbone-only estimate
        try {
            alpha = alphaFromAtmNormalVol<double>(sigma0, problem.forward(), problem.expiry(), o.beta, rho, nu, o.shift);
        } catch (const Error&) {
        }
        z << std::clamp(alpha, lo[0] * 1.001, hi[0] * 0.999), rho, nu;
    }
    return z;
}

// Calibrate one node: bounded Levenberg–Marquardt on S ∘ F with the exact Jacobian, then a verdict.
//
// Bounds. |ρ| < 0.9999 because the expansion has a log singularity at ρ = ±1; ν ≤ 5 because ν²T = 25
// at one year is already far outside the range where a small-time expansion means anything. A fitted
// parameter sitting on either bound is a FAILURE even if the residuals are small: the smile wanted
// something SABR cannot do, and the number is where the solver was stopped, not a calibrated value.
//
// The verdict is taken on RAW residuals (bp of vol), never on the weighted ones, so "max residual
// 0.1 bp" means the same thing on every node whatever the weighting.
SabrFitResult fitSabrSmile(const SabrNodeProblem& problem) {
    const SabrFitOptions& o = problem.options();
    const auto m = static_cast<Eigen::Index>(problem.residualCount());
    const math::Vector S = problem.scaling(math::Vector());
    const auto [lo, hi] = problem.bounds();

    // With α re-solved from the ATM vol, some (ρ, ν) have no admissible α (the ATM equation has no
    // positive root). The solver cannot handle a throw, so such a point returns a large finite
    // residual: the step is rejected and the damping grows.
    const math::ResidualFn residuals = [&](const math::Vector& z) -> math::Vector {
        try {
            const math::Vector r = problem.residuals(z);
            if (r.allFinite()) return r.cwiseProduct(S);
        } catch (const Error&) {
        }
        return math::Vector::Constant(m, kPenaltyBp);
    };
    const math::JacobianFn jacobian = [&](const math::Vector& z) -> math::Matrix {
        try {
            const math::Matrix J = problem.jacobian(z);
            if (J.allFinite()) return S.asDiagonal() * J;
        } catch (const Error&) {
        }
        return math::Matrix::Zero(m, static_cast<Eigen::Index>(problem.parameterCount()));
    };

    math::LmOptions lm = o.lm;
    lm.lowerBounds = lo;
    lm.upperBounds = hi;
    const math::LmResult fit = math::levenbergMarquardt(residuals, jacobian, sabrInitialGuess(problem), lm);

    SabrFitResult res;
    res.status = fit.status;
    res.iterations = fit.iterations;
    auto flag = [&](std::string text) { res.flags.push_back(std::move(text)); };
    auto bp = [](double x) {
        char buf[32];
        std::snprintf(buf, sizeof buf, "%.3f", x);
        return std::string(buf);
    };

    bool evaluated = false;
    try {
        res.params = problem.params(fit.x);
        const math::Vector raw = problem.residuals(fit.x);
        res.residualsBp.assign(raw.data(), raw.data() + raw.size());
        res.maxAbsBp = raw.cwiseAbs().maxCoeff();
        res.rmsBp = std::sqrt(raw.squaredNorm() / static_cast<double>(m));
        evaluated = raw.allFinite();
    } catch (const Error& e) {
        flag(std::string("model cannot be evaluated at the solution: ") + e.message());
    }

    if (evaluated) {
        if (res.maxAbsBp > o.maxResidualBp) flag("max residual " + bp(res.maxAbsBp) + " bp > " + bp(o.maxResidualBp) + " bp");
        for (std::size_t i = 0; i < res.residualsBp.size(); ++i) {
            if (problem.quotes()[i].normalVol + res.residualsBp[i] * 1e-4 <= 0.0) {
                flag("model vol is not positive at " + problem.instrumentIds()[i]);
                break;
            }
        }
        if (!o.alphaFromAtm) {
            const double atmErr = 1e4 * std::fabs(sabrNormalVol<double>(res.params, problem.forward(), problem.forward(), problem.expiry()) - problem.atmNormalVol());
            if (atmErr > o.maxResidualBp) flag("ATM repriced " + bp(atmErr) + " bp off");
        }
    } else if (res.flags.empty()) {
        flag("non-finite residuals at the solution");
    }
    static const char* const names3[] = {"alpha", "rho", "nu"};
    for (Eigen::Index j = 0; j < fit.x.size(); ++j) {
        const char* name = names3[o.alphaFromAtm ? j + 1 : j];
        const double width = hi[j] - lo[j];
        if (fit.x[j] - lo[j] <= kBoundTolerance * width) flag(std::string(name) + " at lower bound");
        if (hi[j] - fit.x[j] <= kBoundTolerance * width) flag(std::string(name) + " at upper bound");
    }
    if (fit.status == math::LmStatus::MaxIterations) flag("solver hit the iteration limit");
    if (fit.status == math::LmStatus::Failed) flag("solver failed");

    res.ok = res.flags.empty();
    if (!res.ok && o.throwOnFailure) {
        std::string msg = "SABR fit failed for node " + problem.id() + ":";
        for (const auto& f : res.flags) msg += " [" + f + "]";
        throw CalibrationError(std::move(msg), res.residualsBp);
    }
    return res;
}

// --- plumbing -------------------------------------------------------------------------------------

std::vector<SmileQuote> smileQuotes(const market::SwaptionVolNodeQuote& node) {
    std::vector<SmileQuote> out;
    out.reserve(node.smile.size() + 1);
    out.push_back({node.forward, node.atmNormalVol});
    for (const auto& p : node.smile) out.push_back({node.forward + p.offset, p.normalVol});
    std::sort(out.begin(), out.end(), [](const SmileQuote& a, const SmileQuote& b) { return a.strike < b.strike; });
    return out;
}

SabrNodeProblem makeNodeProblem(const market::SwaptionVolNodeQuote& node, const SabrFitOptions& options) {
    return SabrNodeProblem(node.id(), smileQuotes(node), node.forward, market::cubeAxisYears(node.expiry), node.atmNormalVol, options);
}

CubeFit fitCube(const market::SwaptionVolQuotes& quotes, const SabrFitOptions& options) {
    std::vector<market::CubeNode<double>> nodes;
    std::vector<std::string> ids;
    std::vector<SabrFitResult> fits;
    for (const auto& e : quotes.expiries()) {
        for (const auto& t : quotes.tenors()) {
            const auto& q = quotes.find(e, t);
            SabrFitResult fit = withContext("fitCube node " + q.id(), [&] { return fitSabrSmile(makeNodeProblem(q, options)); });
            market::CubeNode<double> n;
            n.forward = q.forward;
            n.atmNormalVol = q.atmNormalVol;
            n.sabr = fit.params;
            nodes.push_back(n);
            ids.push_back(q.id());
            fits.push_back(std::move(fit));
        }
    }
    market::SwaptionCube<double> cube(quotes.asOf(), quotes.index(), quotes.expiries(), quotes.tenors(), std::move(nodes), options.beta, options.shift);
    return CubeFit{std::move(cube), std::move(ids), std::move(fits)};
}

} // namespace rl::vol
