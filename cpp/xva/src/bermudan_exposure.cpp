#include "rl/xva/bermudan_exposure.hpp"

#include "rl/core/error.hpp"
#include "rl/math/rng.hpp"

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>

namespace rl::xva {

namespace {

/// Continuation value at one exercise date on a uniform grid in the state.
struct Surface {
    double lo = 0.0, step = 1.0;
    std::vector<double> continuation;

    double operator()(double x) const noexcept {
        const double u = (x - lo) / step;
        const std::size_t n = continuation.size();
        if (u <= 0.0) return continuation.front();
        if (u >= static_cast<double>(n - 1)) return continuation.back();
        const std::size_t i = static_cast<std::size_t>(u);
        const double w = u - static_cast<double>(i);
        return (1.0 - w) * continuation[i] + w * continuation[i + 1];
    }
};

} // namespace

BermudanExposure bermudanExposure(const pricers::BermudanProblem& problem, const BermudanExposureOptions& opt) {
    const std::size_t n = problem.size();
    RL_REQUIRE(n > 0, "bermudanExposure: no exercise dates");
    RL_REQUIRE(opt.paths >= 100, "bermudanExposure: need at least 100 paths");
    RL_REQUIRE(opt.gridPoints >= 51 && opt.gridPoints % 2 == 1 && opt.stdDevs >= 4.0, "bermudanExposure: need an odd grid of at least 51 points over at least 4 standard deviations");
    const std::size_t m = static_cast<std::size_t>(opt.gridPoints);

    // --- Backward induction: C_i on a grid of ±stdDevs standard deviations of the state at date i ---
    std::vector<Surface> surfaces(n);
    std::vector<double> next; // V_{i+1} on date i + 1's grid
    for (std::size_t i = n; i-- > 0;) {
        const double sd = std::sqrt(problem.zeta(i));
        Surface& s = surfaces[i];
        s.lo = -opt.stdDevs * sd;
        s.step = 2.0 * opt.stdDevs * sd / static_cast<double>(m - 1);
        s.continuation.assign(m, 0.0);
        if (i + 1 < n) {
            // Trapezoid rule over date i + 1's grid against the Gaussian transition density. The step
            // is a small fraction of the transition's standard deviation, where the rule is very accurate.
            const Surface& up = surfaces[i + 1];
            const double var = problem.zeta(i + 1) - problem.zeta(i);
            RL_REQUIRE(var > 0.0, "bermudanExposure: the state has no variance between two exercise dates");
            const double norm = up.step / std::sqrt(2.0 * 3.141592653589793 * var);
            for (std::size_t a = 0; a < m; ++a) {
                const double x = s.lo + s.step * static_cast<double>(a);
                double sum = 0.0;
                for (std::size_t b = 0; b < m; ++b) {
                    const double d = up.lo + up.step * static_cast<double>(b) - x;
                    sum += next[b] * std::exp(-0.5 * d * d / var);
                }
                s.continuation[a] = sum * norm;
            }
        }
        next.resize(m);
        for (std::size_t a = 0; a < m; ++a) next[a] = std::max(problem.exerciseValue(i, s.lo + s.step * static_cast<double>(a)), s.continuation[a]);
    }

    BermudanExposure out;
    {
        // Today: one more Gaussian step, from x = 0 to the first exercise date.
        const Surface& first = surfaces[0];
        const double var = problem.zeta(0);
        double sum = 0.0;
        for (std::size_t b = 0; b < m; ++b) {
            const double x = first.lo + first.step * static_cast<double>(b);
            sum += next[b] * std::exp(-0.5 * x * x / var);
        }
        out.price = sum * first.step / std::sqrt(2.0 * 3.141592653589793 * var);
    }

    // --- Forward simulation under the optimal rule ----------------------------------------------------
    out.times.push_back(0.0);
    for (const auto& r : problem.rights()) out.times.push_back(r.time);
    const std::size_t nt = n + 1;
    out.epe.assign(nt, 0.0);
    out.ene.assign(nt, 0.0);
    out.epeStandardError.assign(nt, 0.0);
    out.optionEpe.assign(nt, 0.0);
    out.exercisedBy.assign(nt, 0.0);
    out.epe[0] = out.optionEpe[0] = out.price;

    std::vector<double> sd(n), sumSq(nt, 0.0);
    for (std::size_t i = 0; i < n; ++i) sd[i] = std::sqrt(std::max(problem.zeta(i) - (i ? problem.zeta(i - 1) : 0.0), 0.0));
    double payoff = 0.0, payoffSq = 0.0;
    math::Rng rng(opt.seed);
    for (std::size_t p = 0; p < opt.paths; ++p) {
        double x = 0.0, realised = 0.0;
        bool exercised = false;
        for (std::size_t i = 0; i < n; ++i) {
            x += sd[i] * rng.gaussian();
            const double h = problem.exerciseValue(i, x);
            double value;
            if (exercised) {
                value = h; // the remaining swap
            } else {
                const double hold = surfaces[i](x);
                if (h > 0.0 && h >= hold) {
                    exercised = true;
                    realised = h;
                    value = h;
                } else {
                    value = hold;
                    out.optionEpe[i + 1] += hold;
                }
            }
            const double pos = std::max(value, 0.0);
            out.epe[i + 1] += pos;
            sumSq[i + 1] += pos * pos;
            out.ene[i + 1] += std::min(value, 0.0);
            if (exercised) out.exercisedBy[i + 1] += 1.0;
        }
        payoff += realised;
        payoffSq += realised * realised;
    }
    const double paths = static_cast<double>(opt.paths);
    for (std::size_t i = 1; i < nt; ++i) {
        out.epe[i] /= paths;
        out.ene[i] /= paths;
        out.optionEpe[i] /= paths;
        out.exercisedBy[i] /= paths;
        out.epeStandardError[i] = std::sqrt(std::max(sumSq[i] / paths - out.epe[i] * out.epe[i], 0.0) / paths);
    }
    out.simulatedPrice = payoff / paths;
    out.simulatedPriceStandardError = std::sqrt(std::max(payoffSq / paths - out.simulatedPrice * out.simulatedPrice, 0.0) / paths);
    return out;
}

namespace {

using models::G2pp;

/// What the two-factor simulation needs at one exercise date.
struct Date2F {
    std::vector<double> l1, l2, c;           ///< exercise value: Σ c_k exp(−l1_k x_1 − l2_k x_2)
    double a11 = 0.0, a21 = 0.0, a22 = 0.0;  ///< Cholesky factor of the state increment's covariance
    double sd1 = 1.0, beta = 0.0, sd2 = 0.0; ///< standardisation of the regressors

    double exerciseValue(double x1, double x2) const noexcept {
        double h = 0.0;
        for (std::size_t k = 0; k < c.size(); ++k) h += c[k] * std::exp(-l1[k] * x1 - l2[k] * x2);
        return h;
    }
    /// Monomials u^a v^b, a + b ≤ degree, in the standardised first state and the standardised part of
    /// the second that is independent of it (powers of u alone when there is no such part).
    void basis(int degree, double x1, double x2, std::vector<double>& out) const {
        const double u = x1 / sd1;
        out.clear();
        if (sd2 == 0.0) {
            double pu = 1.0;
            for (int a = 0; a <= degree; ++a, pu *= u) out.push_back(pu);
            return;
        }
        const double v = (x2 - beta * x1) / sd2;
        double pu = 1.0;
        for (int a = 0; a <= degree; ++a, pu *= u) {
            double pv = 1.0;
            for (int b = 0; a + b <= degree; ++b, pv *= v) out.push_back(pu * pv);
        }
    }
};

} // namespace

BermudanExposure bermudanExposureG2pp(const G2pp& model, const std::vector<models::ExerciseRight>& rights, const BermudanExposureG2ppOptions& opt) {
    const std::size_t n = rights.size();
    RL_REQUIRE(n > 0, "bermudanExposureG2pp: no exercise dates");
    for (std::size_t i = 1; i < n; ++i) RL_REQUIRE(rights[i].time > rights[i - 1].time, "bermudanExposureG2pp: rights must be in increasing exercise time");
    RL_REQUIRE(opt.trainingPaths >= 1000 && opt.paths >= 100, "bermudanExposureG2pp: need at least 1000 training paths and 100 valuation paths");
    RL_REQUIRE(opt.degree >= 1 && opt.degree <= 6, "bermudanExposureG2pp: degree must be in [1, 6]");

    // Loadings measured from the last cashflow date (models/g2pp.hpp).
    double horizon = 0.0;
    for (const auto& r : rights)
        for (const auto& f : r.flows) horizon = std::max(horizon, f.time);
    std::vector<Date2F> dates(n);
    G2pp::Covariance prev{};
    for (std::size_t i = 0; i < n; ++i) {
        const G2pp::Covariance cov = model.covariance(rights[i].time);
        Date2F& d = dates[i];
        for (const models::BondFlow& f : rights[i].flows) {
            const G2pp::Loading l = model.loading(f.time, horizon);
            d.l1.push_back(l.h1);
            d.l2.push_back(l.h2);
            d.c.push_back(f.amount * f.df * std::exp(-0.5 * G2pp::quadraticForm(l, cov)));
        }
        const double d11 = std::max(cov.c11 - prev.c11, 0.0), d12 = cov.c12 - prev.c12, d22 = std::max(cov.c22 - prev.c22, 0.0);
        d.a11 = std::sqrt(d11);
        d.a21 = d.a11 > 0.0 ? d12 / d.a11 : 0.0;
        d.a22 = std::sqrt(std::max(d22 - d.a21 * d.a21, 0.0));
        d.sd1 = cov.c11 > 0.0 ? std::sqrt(cov.c11) : 1.0;
        d.beta = cov.c11 > 0.0 ? cov.c12 / cov.c11 : 0.0;
        const double cond = cov.c22 - cov.c12 * d.beta;
        d.sd2 = cond > 1e-14 * std::max(cov.c22, 1e-300) ? std::sqrt(cond) : 0.0;
        prev = cov;
    }

    // --- Fit the continuation value at every date but the last, backwards, on all paths --------------
    const std::size_t N = opt.trainingPaths;
    std::vector<std::vector<double>> x1(n, std::vector<double>(N)), x2(n, std::vector<double>(N)), h(n, std::vector<double>(N));
    {
        math::Rng rng(opt.seed);
        for (std::size_t p = 0; p < N; ++p) {
            double s1 = 0.0, s2 = 0.0;
            for (std::size_t i = 0; i < n; ++i) {
                const double z1 = rng.gaussian(), z2 = rng.gaussian();
                s1 += dates[i].a11 * z1;
                s2 += dates[i].a21 * z1 + dates[i].a22 * z2;
                x1[i][p] = s1;
                x2[i][p] = s2;
                h[i][p] = dates[i].exerciseValue(s1, s2);
            }
        }
    }
    // Two fits per date. In the money, where the decision is made, the fit uses in-the-money paths
    // only: that is the pricing regression, and it puts the exercise boundary in the right place. Out
    // of the money there is nothing to decide but there is still a value, fitted on all paths.
    std::vector<std::vector<double>> inTheMoney(n), everywhere(n);
    std::vector<double> row;
    auto fitted = [&](const std::vector<double>& b, std::size_t i, double s1, double s2) {
        dates[i].basis(opt.degree, s1, s2, row);
        double c = 0.0;
        for (std::size_t k = 0; k < row.size(); ++k) c += b[k] * row[k];
        return c;
    };
    // The value of holding in state (s1, s2) where exercising is worth hv; never negative.
    auto hold = [&](std::size_t i, double s1, double s2, double hv) {
        if (i + 1 == n) return 0.0;
        const bool itm = hv > 0.0 && !inTheMoney[i].empty();
        return std::max(fitted(itm ? inTheMoney[i] : everywhere[i], i, s1, s2), 0.0);
    };
    auto regress = [&](std::size_t i, const std::vector<std::size_t>& rows, const std::vector<double>& target) {
        dates[i].basis(opt.degree, 0.0, 0.0, row);
        const std::size_t terms = row.size();
        Eigen::MatrixXd A(static_cast<Eigen::Index>(rows.size()), static_cast<Eigen::Index>(terms));
        Eigen::VectorXd y(static_cast<Eigen::Index>(rows.size()));
        for (std::size_t r = 0; r < rows.size(); ++r) {
            dates[i].basis(opt.degree, x1[i][rows[r]], x2[i][rows[r]], row);
            for (std::size_t k = 0; k < terms; ++k) A(static_cast<Eigen::Index>(r), static_cast<Eigen::Index>(k)) = row[k];
            y[static_cast<Eigen::Index>(r)] = target[rows[r]];
        }
        const Eigen::VectorXd fit = A.colPivHouseholderQr().solve(y);
        return std::vector<double>(fit.data(), fit.data() + fit.size());
    };
    std::vector<double> cash(N);
    for (std::size_t p = 0; p < N; ++p) cash[p] = std::max(h[n - 1][p], 0.0);
    for (std::size_t i = n - 1; i-- > 0;) {
        std::vector<std::size_t> all(N), itm;
        for (std::size_t p = 0; p < N; ++p) {
            all[p] = p;
            if (h[i][p] > 0.0) itm.push_back(p);
        }
        everywhere[i] = regress(i, all, cash);
        dates[i].basis(opt.degree, 0.0, 0.0, row);
        if (itm.size() >= 10 * row.size()) inTheMoney[i] = regress(i, itm, cash);
        for (std::size_t p : itm)
            if (!inTheMoney[i].empty() && h[i][p] >= hold(i, x1[i][p], x2[i][p], h[i][p])) cash[p] = h[i][p];
    }

    // --- Value the position on fresh paths --------------------------------------------------------------
    BermudanExposure out;
    out.times.push_back(0.0);
    for (const auto& r : rights) out.times.push_back(r.time);
    const std::size_t nt = n + 1;
    out.epe.assign(nt, 0.0);
    out.ene.assign(nt, 0.0);
    out.epeStandardError.assign(nt, 0.0);
    out.optionEpe.assign(nt, 0.0);
    out.exercisedBy.assign(nt, 0.0);
    std::vector<double> sumSq(nt, 0.0);
    double payoff = 0.0, payoffSq = 0.0;
    math::Rng rng(math::SeedPolicy::derive(opt.seed, 1));
    for (std::size_t p = 0; p < opt.paths; ++p) {
        double s1 = 0.0, s2 = 0.0, realised = 0.0;
        bool exercised = false;
        for (std::size_t i = 0; i < n; ++i) {
            const double z1 = rng.gaussian(), z2 = rng.gaussian();
            s1 += dates[i].a11 * z1;
            s2 += dates[i].a21 * z1 + dates[i].a22 * z2;
            const double hv = dates[i].exerciseValue(s1, s2);
            double value;
            if (exercised) {
                value = hv;
            } else {
                const double keep = hold(i, s1, s2, hv);
                if (hv > 0.0 && hv >= keep && (i + 1 == n || !inTheMoney[i].empty())) {
                    exercised = true;
                    realised = hv;
                    value = hv;
                } else {
                    value = keep;
                    out.optionEpe[i + 1] += keep;
                }
            }
            const double pos = std::max(value, 0.0);
            out.epe[i + 1] += pos;
            sumSq[i + 1] += pos * pos;
            out.ene[i + 1] += std::min(value, 0.0);
            if (exercised) out.exercisedBy[i + 1] += 1.0;
        }
        payoff += realised;
        payoffSq += realised * realised;
    }
    const double paths = static_cast<double>(opt.paths);
    for (std::size_t i = 1; i < nt; ++i) {
        out.epe[i] /= paths;
        out.ene[i] /= paths;
        out.optionEpe[i] /= paths;
        out.exercisedBy[i] /= paths;
        out.epeStandardError[i] = std::sqrt(std::max(sumSq[i] / paths - out.epe[i] * out.epe[i], 0.0) / paths);
    }
    out.simulatedPrice = out.price = payoff / paths;
    out.simulatedPriceStandardError = std::sqrt(std::max(payoffSq / paths - out.price * out.price, 0.0) / paths);
    out.epe[0] = out.optionEpe[0] = out.price;
    return out;
}

} // namespace rl::xva
