#include "rl/pricers/bermudan_g2pp.hpp"

#include "rl/core/error.hpp"
#include "rl/math/rng.hpp"
#include "rl/models/g2pp_swaption.hpp"
#include "rl/models/lgm_calibration.hpp"

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>

namespace rl::pricers {

namespace {

using models::G2pp;

/// Everything the simulation needs at one exercise date.
struct DateData {
    // Exercise value: Σ_k c_k exp(−l1_k x_1 − l2_k x_2), loadings measured from the last cashflow date.
    std::vector<double> l1, l2, c;
    // Cholesky factor of the state increment's covariance since the previous date.
    double a11 = 0.0, a21 = 0.0, a22 = 0.0;
    // Standardisation of the regressors: u = x_1/sd1, v = (x_2 − beta x_1)/sd2 (v unused when sd2 = 0).
    double sd1 = 1.0, beta = 0.0, sd2 = 0.0;

    double exerciseValue(double x1, double x2) const noexcept {
        double h = 0.0;
        for (std::size_t k = 0; k < c.size(); ++k) h += c[k] * std::exp(-l1[k] * x1 - l2[k] * x2);
        return h;
    }
};

std::vector<DateData> prepare(const G2pp& model, const std::vector<models::ExerciseRight>& rights) {
    double horizon = 0.0;
    for (const auto& r : rights)
        for (const auto& f : r.flows) horizon = std::max(horizon, f.time);
    std::vector<DateData> out(rights.size());
    G2pp::Covariance prev{};
    for (std::size_t i = 0; i < rights.size(); ++i) {
        const G2pp::Covariance cov = model.covariance(rights[i].time);
        DateData& d = out[i];
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
    return out;
}

/// Monomials u^a v^b with a + b ≤ degree (powers of u alone when the model is one-factor at this date).
void basis(const DateData& d, int degree, double x1, double x2, std::vector<double>& out) {
    const double u = x1 / d.sd1;
    out.clear();
    if (d.sd2 == 0.0) {
        double pu = 1.0;
        for (int a = 0; a <= degree; ++a, pu *= u) out.push_back(pu);
        return;
    }
    const double v = (x2 - d.beta * x1) / d.sd2;
    double pu = 1.0;
    for (int a = 0; a <= degree; ++a, pu *= u) {
        double pv = 1.0;
        for (int b = 0; a + b <= degree; ++b, pv *= v) out.push_back(pu * pv);
    }
}

} // namespace

G2ppLsmcResult priceBermudanLsmcG2pp(const G2pp& model, const std::vector<models::ExerciseRight>& rights, const G2ppLsmcOptions& opt) {
    const std::size_t n = rights.size();
    RL_REQUIRE(n > 0, "priceBermudanLsmcG2pp: no exercise dates");
    for (std::size_t i = 1; i < n; ++i) RL_REQUIRE(rights[i].time > rights[i - 1].time, "priceBermudanLsmcG2pp: rights must be in increasing exercise time");
    RL_REQUIRE(opt.trainingPaths >= 100 && opt.pricingPaths >= 100, "priceBermudanLsmcG2pp: need at least 100 paths in each pass");
    RL_REQUIRE(opt.degree >= 0 && opt.degree <= 6, "priceBermudanLsmcG2pp: degree must be in [0, 6]");
    const std::vector<DateData> dates = prepare(model, rights);

    G2ppLsmcResult res;
    for (const auto& r : rights) res.europeans.push_back(models::g2ppSwaptionPrice(model, r));

    // --- Pass 1: fit the rule ------------------------------------------------------------------
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
    std::vector<std::vector<double>> coefficients(n); // empty = "do not exercise here" (before the last date)
    std::vector<double> row;
    auto continuation = [&](std::size_t i, double s1, double s2) {
        basis(dates[i], opt.degree, s1, s2, row);
        double c = 0.0;
        for (std::size_t k = 0; k < row.size(); ++k) c += coefficients[i][k] * row[k];
        return c;
    };
    auto exercise = [&](std::size_t i, double s1, double s2, double hv) {
        if (hv <= 0.0) return false;
        if (i + 1 == n) return true;
        if (coefficients[i].empty()) return false;
        return hv >= continuation(i, s1, s2);
    };

    std::vector<double> cash(N); // realised payoff of the rule from the current date onwards
    for (std::size_t p = 0; p < N; ++p) cash[p] = std::max(h[n - 1][p], 0.0);
    for (std::size_t i = n - 1; i-- > 0;) {
        std::vector<std::size_t> itm;
        for (std::size_t p = 0; p < N; ++p)
            if (h[i][p] > 0.0) itm.push_back(p);
        basis(dates[i], opt.degree, 0.0, 0.0, row);
        const std::size_t terms = row.size();
        if (itm.size() < 10 * terms) continue; // too few in-the-money paths to fit: never exercise here

        Eigen::MatrixXd A(static_cast<Eigen::Index>(itm.size()), static_cast<Eigen::Index>(terms));
        Eigen::VectorXd y(static_cast<Eigen::Index>(itm.size()));
        for (std::size_t r = 0; r < itm.size(); ++r) {
            const std::size_t p = itm[r];
            basis(dates[i], opt.degree, x1[i][p], x2[i][p], row);
            for (std::size_t k = 0; k < terms; ++k) A(static_cast<Eigen::Index>(r), static_cast<Eigen::Index>(k)) = row[k];
            y[static_cast<Eigen::Index>(r)] = cash[p];
        }
        const Eigen::VectorXd fit = A.colPivHouseholderQr().solve(y);
        coefficients[i].assign(fit.data(), fit.data() + fit.size());

        // Longstaff–Schwartz: the fitted value only DECIDES; the value carried back is the realised payoff.
        for (std::size_t p : itm)
            if (h[i][p] >= continuation(i, x1[i][p], x2[i][p])) cash[p] = h[i][p];
    }
    for (double c : cash) res.inSample += c;
    res.inSample /= static_cast<double>(N);

    // --- Pass 2: value the rule on fresh paths -----------------------------------------------------
    math::Rng rng(math::SeedPolicy::derive(opt.seed, 1));
    std::vector<double> z1(n), z2(n);
    std::vector<std::size_t> exercised(n, 0);
    auto run = [&](double sign) {
        double s1 = 0.0, s2 = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            s1 += sign * dates[i].a11 * z1[i];
            s2 += sign * (dates[i].a21 * z1[i] + dates[i].a22 * z2[i]);
            const double hv = dates[i].exerciseValue(s1, s2);
            if (exercise(i, s1, s2, hv)) {
                ++exercised[i];
                return hv;
            }
        }
        return 0.0;
    };
    const std::size_t samples = opt.antithetic ? opt.pricingPaths / 2 : opt.pricingPaths;
    double sum = 0.0, sumSq = 0.0;
    for (std::size_t p = 0; p < samples; ++p) {
        for (std::size_t i = 0; i < n; ++i) {
            z1[i] = rng.gaussian();
            z2[i] = rng.gaussian();
        }
        const double v = opt.antithetic ? 0.5 * (run(1.0) + run(-1.0)) : run(1.0);
        sum += v;
        sumSq += v * v;
    }
    const double mean = sum / static_cast<double>(samples);
    res.lower = mean;
    res.standardError = std::sqrt(std::max(sumSq / static_cast<double>(samples) - mean * mean, 0.0) / static_cast<double>(samples));
    const double paths = static_cast<double>(opt.antithetic ? 2 * samples : samples);
    for (std::size_t c : exercised) res.exerciseProbability.push_back(static_cast<double>(c) / paths);
    res.coefficients = coefficients;
    res.degree = opt.degree;
    return res;
}

G2ppUpperBoundResult priceBermudanUpperBoundG2pp(const G2pp& model, const std::vector<models::ExerciseRight>& rights, const G2ppLsmcResult& lower,
                                                 const G2ppUpperBoundOptions& opt) {
    const std::size_t n = rights.size();
    RL_REQUIRE(n > 0 && lower.coefficients.size() == n, "priceBermudanUpperBoundG2pp: the exercise rule does not match the rights");
    RL_REQUIRE(opt.outerPaths >= 10 && opt.innerPaths >= 10, "priceBermudanUpperBoundG2pp: need at least 10 outer and inner paths");
    const std::vector<DateData> dates = prepare(model, rights);
    std::vector<double> row;
    auto exercise = [&](std::size_t i, double s1, double s2, double hv) {
        if (hv <= 0.0) return false;
        if (i + 1 == n) return true;
        if (lower.coefficients[i].empty()) return false;
        basis(dates[i], lower.degree, s1, s2, row);
        double c = 0.0;
        for (std::size_t k = 0; k < row.size(); ++k) c += lower.coefficients[i][k] * row[k];
        return hv >= c;
    };
    // Value of following the rule from date `from` on, starting from the state at date from − 1.
    auto follow = [&](std::size_t from, double s1, double s2, math::Rng& rng) {
        for (std::size_t i = from; i < n; ++i) {
            const double z1 = rng.gaussian(), z2 = rng.gaussian();
            s1 += dates[i].a11 * z1;
            s2 += dates[i].a21 * z1 + dates[i].a22 * z2;
            const double hv = dates[i].exerciseValue(s1, s2);
            if (exercise(i, s1, s2, hv)) return hv;
        }
        return 0.0;
    };
    math::Rng outer(opt.seed), inner(math::SeedPolicy::derive(opt.seed, 2));
    double sum = 0.0, sumSq = 0.0;
    for (std::size_t p = 0; p < opt.outerPaths; ++p) {
        double s1 = 0.0, s2 = 0.0, m = 0.0, previousContinuation = 0.0, best = -1e300;
        for (std::size_t i = 0; i < n; ++i) {
            const double z1 = outer.gaussian(), z2 = outer.gaussian();
            s1 += dates[i].a11 * z1;
            s2 += dates[i].a21 * z1 + dates[i].a22 * z2;
            const double h = dates[i].exerciseValue(s1, s2);
            double continuation = 0.0;
            if (i + 1 < n) {
                for (std::size_t q = 0; q < opt.innerPaths; ++q) continuation += follow(i + 1, s1, s2, inner);
                continuation /= static_cast<double>(opt.innerPaths);
            }
            const double L = exercise(i, s1, s2, h) ? h : continuation;
            m += i == 0 ? L - lower.lower : L - previousContinuation;
            previousContinuation = continuation;
            best = std::max(best, h - lower.lower - m);
        }
        best = std::max(best, -lower.lower - m);
        sum += best;
        sumSq += best * best;
    }
    const double paths = static_cast<double>(opt.outerPaths);
    G2ppUpperBoundResult r;
    r.gap = sum / paths;
    r.gapStandardError = std::sqrt(std::max(sumSq / paths - r.gap * r.gap, 0.0) / paths);
    r.upper = lower.lower + r.gap;
    r.standardError = std::sqrt(lower.standardError * lower.standardError + r.gapStandardError * r.gapStandardError);
    return r;
}

CalibratedBermudanG2pp calibrateBermudanG2pp(const models::UsdMarket& market, const instruments::Swap& swap, double firstCallYears, int noticeDays,
                                             const models::G2pp& shape) {
    std::vector<models::ExerciseRight> rights;
    for (core::Date d : models::coterminalExerciseDates(swap, firstCallYears, noticeDays, "USD-SIFMA"))
        rights.push_back(models::makeExerciseRight(swap, market.curves, d));
    models::G2ppStripCalibration cal = models::calibrateG2ppToStrip(shape, models::makeStrip(rights, market.cube));
    return CalibratedBermudanG2pp{std::move(rights), std::move(cal)};
}

} // namespace rl::pricers
