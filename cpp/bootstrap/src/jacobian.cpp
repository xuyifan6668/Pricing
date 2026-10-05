#include "rl/bootstrap/jacobian.hpp"

#include "rl/aad/real.hpp"
#include "rl/aad/tape.hpp"
#include "rl/bootstrap/global_bootstrap.hpp"
#include "rl/bootstrap/residuals.hpp"
#include "rl/core/error.hpp"

#include <cmath>
#include <cstdint>
#include <vector>

namespace rl::bootstrap {

// ∂F/∂z by adjoint differentiation. The residuals are templated on the scalar, so the same code that
// the solver runs on double is recorded once on the tape with the nodes as active inputs; one reverse
// sweep per residual then gives that residual's whole row.
// DESIGN: this replaces a hand-derived chain rule through the interpolation weights (two non-zeros
// per query for linear log-DF, four for monotone convex, plus the collar's active set). The adjoint
// row is exact to rounding for either scheme, costs a small multiple of one residual evaluation, and
// cannot drift out of sync with the interpolator. The sparsity is exact too: a residual that never
// read a node leaves that node's adjoint at exactly zero.
math::Matrix dFdz(const BootstrapProblem& problem, const math::Vector& z) {
    const std::size_t n = problem.parameterCount(), m = problem.residualCount();
    RL_REQUIRE(static_cast<std::size_t>(z.size()) == n, "dFdz: z has the wrong size");
    aad::Tape& tape = aad::Tape::current();
    const aad::Tape::Position start = tape.mark();
    const auto first = static_cast<std::uint32_t>(start.nodes);

    std::vector<aad::Real> zr;
    zr.reserve(n);
    for (std::size_t j = 0; j < n; ++j) zr.push_back(aad::Real::input(z[static_cast<Eigen::Index>(j)]));
    math::Matrix J = math::Matrix::Zero(static_cast<Eigen::Index>(m), static_cast<Eigen::Index>(n));
    try {
        const market::CurveSet<aad::Real> curves = problem.curves(zr);
        std::vector<aad::Real> r;
        r.reserve(m);
        for (const BootstrapInstrument& bi : problem.instruments()) r.push_back(residual(bi, curves, problem.fixings()));
        for (std::size_t i = 0; i < m; ++i) {
            if (!r[i].isActive()) continue;
            tape.clearAdjointsFrom(first);
            tape.seedAdjoint(r[i].node(), 1.0);
            tape.propagate(r[i].node(), first);
            for (std::size_t j = 0; j < n; ++j) J(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j)) = tape.adjoint(zr[j].node());
        }
    } catch (...) {
        tape.rewindTo(start);
        throw;
    }
    tape.rewindTo(start);
    return J;
}

// ∂F/∂q. With rate-form residuals r_i = model_i(z) − q_i, so the matrix is −I; a future is quoted in
// price and its residual is a price difference / 100, hence −1/100. It is diagonal only because every
// quote is a contract term of exactly one instrument.
math::Matrix dFdq(const BootstrapProblem& problem, const math::Vector&) {
    const auto n = static_cast<Eigen::Index>(problem.residualCount());
    math::Matrix D = math::Matrix::Zero(n, n);
    for (Eigen::Index i = 0; i < n; ++i)
        D(i, i) = problem.instruments()[static_cast<std::size_t>(i)].type == market::QuoteType::Future ? -0.01 : -1.0;
    return D;
}

// Zero risk → par risk. At the solution F(z*(q); q) = 0 for every q, so
//     dz/dq = −(∂F/∂z)⁻¹ ∂F/∂q        and        ∂P/∂q = −gᵀ (∂F/∂z)⁻¹ ∂F/∂q,   g = ∂P/∂z.
// Computed adjoint-style: ONE transposed solve (∂F/∂z)ᵀ w = g, then −(∂F/∂q)ᵀ w. Forming dz/dq would
// be n solves; the inverse is never formed. The theorem differentiates the equation, not the
// algorithm, so the result is as accurate as the bootstrap is converged.
math::Vector zeroToParRisk(const math::Vector& dPdz, const math::Matrix& dFdzAtSolution, const math::Matrix& dFdqAtSolution) {
    RL_REQUIRE(dFdzAtSolution.rows() == dFdzAtSolution.cols() && dFdzAtSolution.rows() == dPdz.size(), "zeroToParRisk: dimension mismatch");
    Eigen::PartialPivLU<math::Matrix> lu(dFdzAtSolution.transpose());
    if (!(lu.rcond() > 1e-14)) throw NumericalError("zeroToParRisk: the bootstrap Jacobian is singular");
    const math::Vector w = lu.solve(dPdz);
    return -(dFdqAtSolution.transpose() * w);
}

// --- helpers ------------------------------------------------------------------------------------------

math::Vector zeroRiskByBump(const BootstrapProblem& problem, const math::Vector& z, const PvFunction& pv, double step) {
    RL_REQUIRE(step > 0.0, "zeroRiskByBump: step must be positive");
    const auto n = static_cast<Eigen::Index>(problem.parameterCount());
    RL_REQUIRE(z.size() == n, "zeroRiskByBump: z has the wrong size");
    math::Vector g(n);
    math::Vector zp = z, zm = z;
    for (Eigen::Index j = 0; j < n; ++j) {
        const double h = step * (1.0 + std::fabs(z[j]));
        zp[j] = z[j] + h;
        zm[j] = z[j] - h;
        g[j] = (pv(problem.curves(zp)) - pv(problem.curves(zm))) / (2.0 * h);
        zp[j] = z[j];
        zm[j] = z[j];
    }
    return g;
}

math::Vector parRisk(const BootstrapProblem& problem, const math::Vector& z, const math::Vector& dPdz) {
    return zeroToParRisk(dPdz, dFdz(problem, z), dFdq(problem, z));
}

math::Vector parRiskByRebuild(const BootstrapProblem& problem, const math::Vector& z, const PvFunction& pv, double bump) {
    RL_REQUIRE(bump != 0.0, "parRiskByRebuild: bump must be non-zero");
    const auto n = static_cast<Eigen::Index>(problem.residualCount());
    const double base = pv(problem.curves(z));
    const math::Vector q = problem.quoteValues();
    math::Vector out(n);
    BootstrapOptions opt;
    opt.toleranceBp = 1e-10;
    for (Eigen::Index i = 0; i < n; ++i) {
        math::Vector qb = q;
        qb[i] += bump;
        BootstrapProblem bumped = problem.withQuotes(qb);
        BootstrapResult r = solve(bumped, opt);
        RL_REQUIRE(r.converged, "parRiskByRebuild: bumped bootstrap did not converge for " + problem.instrumentIds()[static_cast<std::size_t>(i)]);
        out[i] = (pv(bumped.curves(r.z)) - base) / bump;
    }
    return out;
}

} // namespace rl::bootstrap
