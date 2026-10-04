#include "rl/bootstrap/global_bootstrap.hpp"

#include "rl/bootstrap/jacobian.hpp"
#include "rl/bootstrap/residuals.hpp"
#include "rl/core/error.hpp"

#include <algorithm>
#include <cmath>
#include <string>
#include <string_view>
#include <vector>

namespace rl::bootstrap {

namespace {

/// The quoted rate of an instrument as a decimal rate, for seeding a flat curve.
double quotedRate(const BootstrapInstrument& bi) { return bi.type == market::QuoteType::Future ? 1.0 - bi.quote / 100.0 : bi.quote; }

/// Rate of the shortest outright (non-basis, non-turn) instrument on a curve; nullopt-like −1 flag via bool.
bool frontRate(const BootstrapProblem& p, std::string_view curveName, double& rate) {
    for (const BootstrapInstrument* bi : instrumentsOn(p.instruments(), curveName)) {
        if (bi->type == market::QuoteType::Turn || bi->type == market::QuoteType::BasisSwap) continue;
        rate = quotedRate(*bi);
        return true;
    }
    return false;
}

std::string parameterLabel(const BootstrapProblem& p, std::size_t j) {
    std::size_t off = 0;
    for (const CurveSpec& s : p.curveSpecs()) {
        if (j < off + s.pillarIds.size()) return s.name + " node " + s.pillarIds[j - off];
        off += s.pillarIds.size();
        if (j < off + s.turns.size()) return s.name + " turn " + s.turns[j - off].id;
        off += s.turns.size();
    }
    return "parameter " + std::to_string(j);
}

std::vector<std::size_t> allIndices(std::size_t n) {
    std::vector<std::size_t> v(n);
    for (std::size_t i = 0; i < n; ++i) v[i] = i;
    return v;
}

math::Vector pick(const math::Vector& v, const std::vector<std::size_t>& idx) {
    math::Vector out(static_cast<Eigen::Index>(idx.size()));
    for (std::size_t k = 0; k < idx.size(); ++k) out[static_cast<Eigen::Index>(k)] = v[static_cast<Eigen::Index>(idx[k])];
    return out;
}

/// Damped Newton on the square subsystem {residual rows} = 0 over {parameter cols}, all other
/// parameters frozen at their values in z. The global solve passes every row and column; the
/// sequential bootstrap passes one curve's block at a time.
BootstrapResult newton(const BootstrapProblem& problem, math::Vector z, const std::vector<std::size_t>& rows, const std::vector<std::size_t>& cols,
                       const BootstrapOptions& opt) {
    RL_REQUIRE(rows.size() == cols.size(), "bootstrap: the Newton block must be square");
    const auto m = static_cast<Eigen::Index>(rows.size());
    BootstrapResult res;
    const math::Vector frozenScale = problem.scaling(z);
    auto scaled = [&](const math::Vector& at) -> math::Vector {
        const math::Vector s = opt.rescaleEachIteration ? problem.scaling(at) : frozenScale;
        return pick(problem.residuals(at).cwiseProduct(s), rows);
    };

    math::Vector sf = scaled(z);
    for (int it = 0;; ++it) {
        const double worst = sf.cwiseAbs().maxCoeff();
        if (opt.verbose) res.log.push_back("iter " + std::to_string(it) + "  max |S F| = " + std::to_string(worst) + " bp");
        if (worst < opt.toleranceBp) {
            res.converged = true;
            res.iterations = it;
            break;
        }
        if (it == opt.maxIterations) {
            res.iterations = it;
            res.log.push_back("not converged after " + std::to_string(it) + " iterations; worst residual " + std::to_string(worst) + " bp");
            break;
        }

        // Scaled Jacobian block. S is evaluated at z; with rate-form residuals it is constant.
        const math::Matrix full = opt.analyticJacobian ? dFdz(problem, z) : finiteDifferenceJacobian(problem, z, opt.fdStep);
        const math::Vector s = opt.rescaleEachIteration ? problem.scaling(z) : frozenScale;
        math::Matrix J(m, m);
        for (Eigen::Index i = 0; i < m; ++i)
            for (Eigen::Index j = 0; j < m; ++j)
                J(i, j) = s[static_cast<Eigen::Index>(rows[static_cast<std::size_t>(i)])] *
                          full(static_cast<Eigen::Index>(rows[static_cast<std::size_t>(i)]), static_cast<Eigen::Index>(cols[static_cast<std::size_t>(j)]));

        // DESIGN: LU with partial pivoting, never an explicit inverse. A (near-)singular J means some
        // node is not pinned by any instrument; the column with the smallest norm names it.
        Eigen::PartialPivLU<math::Matrix> lu(J);
        if (!(lu.rcond() > 1e-14)) {
            Eigen::Index jmin = 0;
            J.colwise().norm().minCoeff(&jmin);
            throw NumericalError("bootstrap: singular Jacobian; no instrument pins " + parameterLabel(problem, cols[static_cast<std::size_t>(jmin)]));
        }
        math::Vector step = lu.solve(-sf);
        // Clamp: from a poor guess the linearisation overshoots; no node moves more than maxStep in y.
        for (Eigen::Index j = 0; j < m; ++j) step[j] = std::clamp(step[j], -opt.maxStep, opt.maxStep);

        // Backtracking line search on ‖S F‖₂.
        const double norm0 = sf.norm();
        double lambda = 1.0;
        math::Vector trial = z, sfTrial;
        bool accepted = false;
        for (int h = 0; h <= 20; ++h, lambda *= 0.5) {
            trial = z;
            for (Eigen::Index j = 0; j < m; ++j) trial[static_cast<Eigen::Index>(cols[static_cast<std::size_t>(j)])] += lambda * step[j];
            sfTrial = scaled(trial);
            if (!opt.lineSearch || (sfTrial.allFinite() && sfTrial.norm() < norm0)) {
                accepted = true;
                break;
            }
        }
        if (!accepted) {
            // The residual is at the rounding floor of the par-rate arithmetic (or the quote set is
            // inconsistent): no step along the Newton direction improves it. Report, do not loop.
            Eigen::Index imax = 0;
            sf.cwiseAbs().maxCoeff(&imax);
            res.iterations = it;
            res.log.push_back("stalled: no decrease along the Newton direction; worst residual " + std::to_string(worst) + " bp on " +
                              problem.instrumentIds()[rows[static_cast<std::size_t>(imax)]]);
            break;
        }
        z = trial;
        sf = sfTrial;
    }

    res.z = z;
    res.residuals = problem.residuals(z);
    res.scaledResiduals = res.residuals.cwiseProduct(problem.scaling(z));
    res.maxAbsBp = res.scaledResiduals.cwiseAbs().maxCoeff();
    return res;
}

} // namespace

// Flat curves at the front rate of each curve: y_i = r₀ t_i. A curve with no outright instrument at
// the front (Fed Funds here, pinned by basis swaps) starts from the reference curve's front rate plus
// its shortest spread. Zero would make every discount factor 1 and the first Newton step enormous;
// the par rates themselves are not zero rates and are a worse guess on a steep curve. Production
// would seed from yesterday's solution, which has the same layout: pass it as z.
math::Vector initialGuess(const BootstrapProblem& problem) {
    math::Vector z(static_cast<Eigen::Index>(problem.parameterCount()));
    Eigen::Index k = 0;
    for (const CurveSpec& spec : problem.curveSpecs()) {
        double r0 = 0.03;
        if (!frontRate(problem, spec.name, r0)) {
            for (const BootstrapInstrument* bi : instrumentsOn(problem.instruments(), spec.name)) {
                if (bi->type != market::QuoteType::BasisSwap) continue;
                double ref = 0.03;
                frontRate(problem, bi->source.referenceIndex, ref);
                r0 = ref + bi->quote;
                break;
            }
        }
        for (double t : spec.pillarTimes) z[k++] = r0 * t;
        for (const TurnOfYear& turn : spec.turns) z[k++] = turn.jump;
    }
    return z;
}

BootstrapResult solve(const BootstrapProblem& problem, const BootstrapOptions& options) {
    const std::vector<std::size_t> all = allIndices(problem.parameterCount());
    BootstrapResult r = newton(problem, initialGuess(problem), all, all, options);
    return r;
}

// One curve at a time, in curveSpecs() order, each solved against whatever the other curves currently
// are. Kept to show what goes wrong: the 25Y/30Y basis swaps that pin the SOFR long end are solved
// against the *guessed* Fed Funds curve, and are left mispriced once Fed Funds is solved afterwards.
BootstrapResult solveSequential(const BootstrapProblem& problem, const BootstrapOptions& options) {
    math::Vector z = initialGuess(problem);
    int iterations = 0;
    std::vector<std::string> log;
    for (const CurveSpec& spec : problem.curveSpecs()) {
        const auto [off, count] = problem.parameterRange(spec.name);
        std::vector<std::size_t> cols(count), rows;
        for (std::size_t j = 0; j < count; ++j) cols[j] = off + j;
        for (std::size_t i = 0; i < problem.size(); ++i)
            if (problem.instruments()[i].curveName == spec.name) rows.push_back(i);
        BootstrapResult block = newton(problem, z, rows, cols, options);
        z = block.z;
        iterations += block.iterations;
        log.insert(log.end(), block.log.begin(), block.log.end());
    }
    BootstrapResult res;
    res.z = z;
    res.residuals = problem.residuals(z);
    res.scaledResiduals = res.residuals.cwiseProduct(problem.scaling(z));
    res.maxAbsBp = res.scaledResiduals.cwiseAbs().maxCoeff();
    res.iterations = iterations;
    res.converged = res.maxAbsBp < options.toleranceBp;
    res.log = std::move(log);
    return res;
}

// --- helpers ------------------------------------------------------------------------------------------

math::Matrix finiteDifferenceJacobian(const BootstrapProblem& problem, const math::Vector& z, double step) {
    RL_REQUIRE(step > 0.0, "finiteDifferenceJacobian: step must be positive");
    const auto n = static_cast<Eigen::Index>(problem.parameterCount());
    const auto m = static_cast<Eigen::Index>(problem.residualCount());
    RL_REQUIRE(z.size() == n, "finiteDifferenceJacobian: z has the wrong size");
    math::Matrix J(m, n);
    math::Vector zp = z, zm = z;
    for (Eigen::Index j = 0; j < n; ++j) {
        const double h = step * (1.0 + std::fabs(z[j]));
        zp[j] = z[j] + h;
        zm[j] = z[j] - h;
        J.col(j) = (problem.residuals(zp) - problem.residuals(zm)) / (2.0 * h);
        zp[j] = z[j];
        zm[j] = z[j];
    }
    return J;
}

math::Vector scaledResiduals(const BootstrapProblem& problem, const math::Vector& z) {
    const math::Vector r = problem.residuals(z);
    const math::Vector s = problem.scaling(z);
    RL_REQUIRE(r.size() == s.size(), "scaledResiduals: residual/scale size mismatch");
    return r.cwiseProduct(s);
}

} // namespace rl::bootstrap
