// Gauss–Hermite nodes and weights by Golub–Welsch (1969) — eigenvalues of the
// symmetric tridiagonal Jacobi matrix with off-diagonal √(k/2) — normalised to E[g(Z)] for Z ~ N(0,1).
#include "rl/math/quadrature.hpp"

#include "rl/core/error.hpp"

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <numeric>

namespace rl::math {

GaussHermite::GaussHermite(std::size_t n) {
    RL_REQUIRE(n >= 1 && n <= 512, "Gauss-Hermite order must be in [1, 512]");
    build(n);
}

void GaussHermite::build(std::size_t n) {
    x_.assign(n, 0.0);
    w_.assign(n, 0.0);
    if (n == 1) {
        w_[0] = 1.0;
        return;
    }
    const auto m = static_cast<Eigen::Index>(n);
    Eigen::MatrixXd J = Eigen::MatrixXd::Zero(m, m);
    for (Eigen::Index k = 1; k < m; ++k) J(k, k - 1) = J(k - 1, k) = std::sqrt(0.5 * static_cast<double>(k));
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(J);
    RL_REQUIRE(es.info() == Eigen::Success, "Gauss-Hermite: eigen-decomposition failed");
    // Physicists' rule: ξ_i = eigenvalues, ω_i = √π · v_{0i}²; probabilists' form x = √2 ξ, w = ω/√π.
    std::vector<std::size_t> order(n);
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) { return es.eigenvalues()[static_cast<Eigen::Index>(a)] < es.eigenvalues()[static_cast<Eigen::Index>(b)]; });
    for (std::size_t i = 0; i < n; ++i) {
        const auto k = static_cast<Eigen::Index>(order[i]);
        x_[i] = std::sqrt(2.0) * es.eigenvalues()[k];
        const double v0 = es.eigenvectors()(0, k);
        w_[i] = v0 * v0;
    }
    // Symmetrise to kill 1e-16 asymmetry and renormalise the weights to sum exactly to one.
    for (std::size_t i = 0; i < n / 2; ++i) {
        const double xa = 0.5 * (x_[n - 1 - i] - x_[i]);
        x_[i] = -xa;
        x_[n - 1 - i] = xa;
        const double wa = 0.5 * (w_[i] + w_[n - 1 - i]);
        w_[i] = w_[n - 1 - i] = wa;
    }
    if (n % 2 == 1) x_[n / 2] = 0.0;
    const double sum = std::accumulate(w_.begin(), w_.end(), 0.0);
    for (double& w : w_) w /= sum;
}

} // namespace rl::math
