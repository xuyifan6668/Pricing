// Levenberg–Marquardt with Marquardt scaling and box clamping (Nocedal &
// Wright, Numerical Optimization 2e, §10.3; Marquardt 1963), plus the central-difference Jacobian.
#include "rl/math/lm.hpp"

#include "rl/core/error.hpp"

#include <cmath>

namespace rl::math {

const char* toString(LmStatus s) {
    switch (s) {
    case LmStatus::Converged: return "Converged";
    case LmStatus::GradientSmall: return "GradientSmall";
    case LmStatus::StepSmall: return "StepSmall";
    case LmStatus::CostSmall: return "CostSmall";
    case LmStatus::MaxIterations: return "MaxIterations";
    case LmStatus::Failed: return "Failed";
    }
    return "?";
}

Matrix finiteDifferenceJacobian(const ResidualFn& residuals, const Vector& x, double relStep) {
    RL_REQUIRE(relStep > 0.0, "finiteDifferenceJacobian: relStep must be positive");
    const Eigen::Index n = x.size();
    Vector xp = x, xm = x;
    Matrix J;
    for (Eigen::Index j = 0; j < n; ++j) {
        const double h = relStep * (1.0 + std::fabs(x[j]));
        xp[j] = x[j] + h;
        xm[j] = x[j] - h;
        const Vector rp = residuals(xp), rm = residuals(xm);
        RL_REQUIRE(rp.size() == rm.size(), "finiteDifferenceJacobian: residual size changed");
        if (J.size() == 0) J.resize(rp.size(), n);
        J.col(j) = (rp - rm) / (2.0 * h);
        xp[j] = x[j];
        xm[j] = x[j];
    }
    if (J.size() == 0) J.resize(residuals(x).size(), 0);
    return J;
}

namespace {
Vector clampToBox(Vector x, const LmOptions& opt) {
    if (opt.lowerBounds.size() == x.size()) x = x.cwiseMax(opt.lowerBounds);
    if (opt.upperBounds.size() == x.size()) x = x.cwiseMin(opt.upperBounds);
    return x;
}
} // namespace

LmResult levenbergMarquardt(const ResidualFn& residuals, const JacobianFn& jac, const Vector& x0, LmOptions opt) {
    RL_REQUIRE(opt.maxIter >= 0 && opt.lambda0 > 0.0, "levenbergMarquardt: bad options");
    auto jacobian = [&](const Vector& x) { return jac ? jac(x) : finiteDifferenceJacobian(residuals, x); };

    LmResult res;
    res.x = clampToBox(x0, opt);
    res.residuals = residuals(res.x);
    ++res.evaluations;
    res.cost = 0.5 * res.residuals.squaredNorm();
    const Eigen::Index n = res.x.size();
    double lambda = opt.lambda0;
    Matrix J = jacobian(res.x);
    RL_REQUIRE(J.rows() == res.residuals.size() && J.cols() == n, "levenbergMarquardt: Jacobian shape mismatch");

    for (int iter = 0;; ++iter) {
        const Vector g = J.transpose() * res.residuals;
        if (g.size() == 0 || g.lpNorm<Eigen::Infinity>() <= opt.gradTol) {
            res.status = LmStatus::GradientSmall;
            return res;
        }
        if (iter >= opt.maxIter) {
            res.status = LmStatus::MaxIterations;
            return res;
        }
        ++res.iterations;
        const Matrix JtJ = J.transpose() * J;
        Vector D = opt.useMarquardtScaling ? Vector(JtJ.diagonal()) : Vector::Ones(n);
        const double dmax = D.maxCoeff();
        D = D.cwiseMax(1e-12 * (dmax > 0.0 ? dmax : 1.0));

        bool accepted = false;
        for (int attempt = 0; attempt < 60 && !accepted; ++attempt) {
            const Matrix A = JtJ + lambda * D.asDiagonal().toDenseMatrix();
            Eigen::LDLT<Matrix> ldlt(A);
            if (ldlt.info() != Eigen::Success) {
                lambda *= opt.lambdaUp;
                continue;
            }
            const Vector step = ldlt.solve(-g);
            const Vector xNew = clampToBox(res.x + step, opt);
            const Vector rNew = residuals(xNew);
            ++res.evaluations;
            const double costNew = 0.5 * rNew.squaredNorm();
            if (costNew < res.cost) {
                const double stepNorm = (xNew - res.x).norm();
                const double costDrop = res.cost - costNew;
                res.x = xNew;
                res.residuals = rNew;
                const double oldCost = res.cost;
                res.cost = costNew;
                lambda = std::max(lambda * opt.lambdaDown, 1e-300);
                accepted = true;
                if (stepNorm <= opt.stepTol * (1.0 + res.x.norm())) {
                    res.status = LmStatus::StepSmall;
                    return res;
                }
                if (oldCost > 0.0 && costDrop <= opt.costTol * oldCost) {
                    res.status = LmStatus::CostSmall;
                    return res;
                }
                if (res.cost <= 1e-300) {
                    res.status = LmStatus::Converged;
                    return res;
                }
            } else {
                lambda *= opt.lambdaUp;
            }
        }
        if (!accepted) {
            res.status = LmStatus::StepSmall; // damping has grown until the step vanished: a local minimum
            return res;
        }
        J = jacobian(res.x);
    }
}

} // namespace rl::math
