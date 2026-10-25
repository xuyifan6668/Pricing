#include "rl/pricers/bermudan_mc.hpp"

#include "rl/core/error.hpp"
#include "rl/math/rng.hpp"

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>

namespace rl::pricers {

namespace {

/// Standard deviations of the state increments between exercise dates.
std::vector<double> stepDeviations(const BermudanProblem& p) {
    std::vector<double> sd(p.size());
    double prev = 0.0;
    for (std::size_t i = 0; i < p.size(); ++i) {
        const double z = p.zeta(i);
        sd[i] = std::sqrt(std::max(z - prev, 0.0));
        prev = z;
    }
    return sd;
}

double regressor(const ExercisePolicy& policy, std::size_t i, double x, double h) {
    return (policy.basis == RegressionBasis::State ? x : h) / policy.scale[i];
}

/// Value of following the policy from date `from` onwards, starting from state x just before
/// date `from` is simulated (i.e. x is the state at date from − 1, or 0 at time 0).
double followPolicy(const BermudanProblem& p, const ExercisePolicy& policy, const std::vector<double>& sd, std::size_t from, double x, math::Rng& rng) {
    for (std::size_t i = from; i < p.size(); ++i) {
        x += sd[i] * rng.gaussian();
        const double h = p.exerciseValue(i, x);
        if (policy.exercise(p, i, x, h)) return h;
    }
    return 0.0;
}

} // namespace

double ExercisePolicy::continuation(const BermudanProblem&, std::size_t i, double x, double h) const {
    const std::vector<double>& b = coefficients[i];
    const double f = regressor(*this, i, x, h);
    double c = 0.0;
    for (std::size_t k = b.size(); k-- > 0;) c = c * f + b[k]; // Horner
    return c;
}

bool ExercisePolicy::exercise(const BermudanProblem& p, std::size_t i, double x, double h) const {
    if (h <= 0.0) return false;
    if (i + 1 == p.size()) return true;
    if (coefficients[i].empty()) return false;
    return h >= continuation(p, i, x, h);
}

LsmcResult priceBermudanLsmc(const BermudanProblem& problem, const LsmcOptions& opt) {
    const std::size_t n = problem.size();
    RL_REQUIRE(n > 0, "priceBermudanLsmc: no exercise dates");
    RL_REQUIRE(opt.trainingPaths >= 100 && opt.pricingPaths >= 100, "priceBermudanLsmc: need at least 100 paths in each pass");
    RL_REQUIRE(opt.degree >= 0 && opt.degree <= 8, "priceBermudanLsmc: degree must be in [0, 8]");
    const std::vector<double> sd = stepDeviations(problem);
    const std::size_t terms = static_cast<std::size_t>(opt.degree) + 1;

    LsmcResult res;
    res.policy.basis = opt.basis;
    res.policy.degree = opt.degree;
    res.policy.coefficients.assign(n, {});
    res.policy.scale.assign(n, 1.0);

    // --- Pass 1: fit the rule ------------------------------------------------------------------
    const std::size_t N = opt.trainingPaths;
    std::vector<std::vector<double>> x(n, std::vector<double>(N)), h(n, std::vector<double>(N));
    {
        math::Rng rng(opt.seed);
        for (std::size_t p = 0; p < N; ++p) {
            double s = 0.0;
            for (std::size_t i = 0; i < n; ++i) {
                s += sd[i] * rng.gaussian();
                x[i][p] = s;
                h[i][p] = problem.exerciseValue(i, s);
            }
        }
    }
    std::vector<double> cash(N); // realised payoff of the rule from the current date onwards
    for (std::size_t p = 0; p < N; ++p) cash[p] = std::max(h[n - 1][p], 0.0);
    for (std::size_t i = n - 1; i-- > 0;) {
        std::vector<std::size_t> itm;
        for (std::size_t p = 0; p < N; ++p)
            if (h[i][p] > 0.0) itm.push_back(p);
        if (itm.size() < 10 * terms) continue; // too few in-the-money paths to fit: never exercise here

        // Normalise the regressor so the powers are O(1) and the least-squares problem is well conditioned.
        double scale = 0.0;
        for (std::size_t p : itm) scale += std::fabs(opt.basis == RegressionBasis::State ? x[i][p] : h[i][p]);
        scale /= static_cast<double>(itm.size());
        res.policy.scale[i] = scale > 0.0 ? scale : 1.0;

        Eigen::MatrixXd A(static_cast<Eigen::Index>(itm.size()), static_cast<Eigen::Index>(terms));
        Eigen::VectorXd y(static_cast<Eigen::Index>(itm.size()));
        for (std::size_t r = 0; r < itm.size(); ++r) {
            const std::size_t p = itm[r];
            const double f = regressor(res.policy, i, x[i][p], h[i][p]);
            double pw = 1.0;
            for (std::size_t k = 0; k < terms; ++k, pw *= f) A(static_cast<Eigen::Index>(r), static_cast<Eigen::Index>(k)) = pw;
            y[static_cast<Eigen::Index>(r)] = cash[p];
        }
        const Eigen::VectorXd beta = A.colPivHouseholderQr().solve(y);
        res.policy.coefficients[i].assign(beta.data(), beta.data() + beta.size());

        // Longstaff–Schwartz: the fitted value only DECIDES; the value carried back is the realised payoff.
        for (std::size_t p : itm)
            if (h[i][p] >= res.policy.continuation(problem, i, x[i][p], h[i][p])) cash[p] = h[i][p];
    }
    for (double c : cash) res.inSample += c;
    res.inSample /= static_cast<double>(N);

    // --- Pass 2: value the rule on fresh paths -----------------------------------------------------
    // DESIGN: antithetic pairs only. A control variate on the best European (whose price is known in
    // closed form) was tried and cut the standard error by just 12% on top of antithetic sampling,
    // which already removes the part of the payoff that is odd in the path; not worth the code.
    math::Rng rng(math::SeedPolicy::derive(opt.seed, 1));
    std::vector<double> z(n);
    std::vector<std::size_t> exercised(n, 0);
    auto run = [&](double sign) {
        double s = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            s += sign * sd[i] * z[i];
            const double hv = problem.exerciseValue(i, s);
            if (res.policy.exercise(problem, i, s, hv)) {
                ++exercised[i];
                return hv;
            }
        }
        return 0.0;
    };
    const std::size_t samples = opt.antithetic ? opt.pricingPaths / 2 : opt.pricingPaths;
    double sum = 0.0, sumSq = 0.0;
    for (std::size_t p = 0; p < samples; ++p) {
        for (std::size_t i = 0; i < n; ++i) z[i] = rng.gaussian();
        const double v = opt.antithetic ? 0.5 * (run(1.0) + run(-1.0)) : run(1.0);
        sum += v;
        sumSq += v * v;
    }
    const double mean = sum / static_cast<double>(samples);
    res.lower = mean;
    res.standardError = std::sqrt(std::max(sumSq / static_cast<double>(samples) - mean * mean, 0.0) / static_cast<double>(samples));
    const double paths = static_cast<double>(opt.antithetic ? 2 * samples : samples);
    for (std::size_t c : exercised) res.exerciseProbability.push_back(static_cast<double>(c) / paths);
    return res;
}

UpperBoundResult priceBermudanUpperBound(const BermudanProblem& problem, const LsmcResult& lower, const UpperBoundOptions& opt) {
    const ExercisePolicy& policy = lower.policy;
    const std::size_t n = problem.size();
    RL_REQUIRE(n > 0 && policy.coefficients.size() == n, "priceBermudanUpperBound: policy does not match the problem");
    RL_REQUIRE(opt.outerPaths >= 10 && opt.innerPaths >= 10, "priceBermudanUpperBound: need at least 10 outer and inner paths");
    const std::vector<double> sd = stepDeviations(problem);
    math::Rng outer(opt.seed), inner(math::SeedPolicy::derive(opt.seed, 2));

    double sum = 0.0, sumSq = 0.0;
    for (std::size_t p = 0; p < opt.outerPaths; ++p) {
        // M is the martingale part of L, the value of following the rule: M_0 = L_0 and
        //     M_i = M_{i−1} + L_i − E_{i−1}[L_i].
        // m below is M_i − L_0, so the pathwise dual value max_i (h_i − m_i) is the gap contribution.
        // "Never exercise" (payoff 0) is one of the holder's choices at the end.
        double x = 0.0, m = 0.0, previousContinuation = 0.0, best = -1e300;
        for (std::size_t i = 0; i < n; ++i) {
            x += sd[i] * outer.gaussian();
            const double h = problem.exerciseValue(i, x);
            // E_i[value of following the rule from date i + 1]: nested simulation from this state.
            double continuation = 0.0;
            if (i + 1 < n) {
                for (std::size_t q = 0; q < opt.innerPaths; ++q) continuation += followPolicy(problem, policy, sd, i + 1, x, inner);
                continuation /= static_cast<double>(opt.innerPaths);
            }
            const double L = policy.exercise(problem, i, x, h) ? h : continuation;
            // First date: E_0[L_1] = L_0, which is what m is measured from.
            m += i == 0 ? L - lower.lower : L - previousContinuation;
            previousContinuation = continuation;
            best = std::max(best, h - lower.lower - m);
        }
        best = std::max(best, -lower.lower - m);
        sum += best;
        sumSq += best * best;
    }
    const double paths = static_cast<double>(opt.outerPaths);
    UpperBoundResult res;
    res.gap = sum / paths;
    res.gapStandardError = std::sqrt(std::max(sumSq / paths - res.gap * res.gap, 0.0) / paths);
    res.upper = lower.lower + res.gap;
    res.standardError = std::sqrt(lower.standardError * lower.standardError + res.gapStandardError * res.gapStandardError);
    return res;
}

} // namespace rl::pricers
