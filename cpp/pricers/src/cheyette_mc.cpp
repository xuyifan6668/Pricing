#include "rl/pricers/cheyette_mc.hpp"

#include "rl/core/error.hpp"
#include "rl/math/rng.hpp"
#include "rl/math/rootfind.hpp"
#include "rl/models/lgm_swaption.hpp"

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>

namespace rl::pricers {

namespace {

using models::Cheyette;

/// The state of every path: x, y and the running integral of x. Path p and path p + n/2 are an
/// antithetic pair.
struct State {
    std::vector<double> x, y, integral;
    explicit State(std::size_t n) : x(n, 0.0), y(n, 0.0), integral(n, 0.0) {}
    std::size_t size() const noexcept { return x.size(); }
};

int stepsBetween(double t0, double t1, int stepsPerYear) { return std::max(1, static_cast<int>(std::ceil((t1 - t0) * stepsPerYear - 1e-9))); }

/// Standard normals for one interval: `steps` rows of n/2 draws (the other half is their negative).
std::vector<double> drawNormals(math::Rng& rng, std::size_t half, int steps) {
    std::vector<double> z(half * static_cast<std::size_t>(steps));
    for (double& v : z) v = rng.gaussian();
    return z;
}

/// Euler steps from t0 to t1 with volatility level `level` (constant on the interval).
void advance(State& s, const Cheyette& model, double level, double t0, double t1, int steps, const std::vector<double>& normals) {
    const std::size_t n = s.size(), half = n / 2;
    const double dt = (t1 - t0) / steps, root = std::sqrt(dt), kappa = model.kappa();
    for (int k = 0; k < steps; ++k) {
        const double* z = normals.data() + static_cast<std::size_t>(k) * half;
        for (std::size_t p = 0; p < n; ++p) {
            const double dw = (p < half ? z[p] : -z[p - half]) * root;
            const double x0 = s.x[p], eta = level * model.shape(x0);
            const double x1 = x0 + (s.y[p] - kappa * x0) * dt + eta * dw;
            s.y[p] += (eta * eta - 2.0 * kappa * s.y[p]) * dt;
            s.integral[p] += 0.5 * (x0 + x1) * dt;
            s.x[p] = x1;
        }
    }
}

/// A right's discounted exercise value as a function of the state at its exercise time.
struct Payoff {
    std::vector<double> g, c; ///< G(T_k − t) and a_k D(T_k)
    Payoff(const Cheyette& model, const models::ExerciseRight& r) {
        for (const models::BondFlow& f : r.flows) {
            g.push_back(model.G(f.time - r.time));
            c.push_back(f.amount * f.df);
        }
    }
    double operator()(double x, double y, double integral) const noexcept {
        double v = 0.0;
        for (std::size_t k = 0; k < c.size(); ++k) v += c[k] * std::exp(-g[k] * x - 0.5 * g[k] * g[k] * y);
        return std::exp(-integral) * v;
    }
};

/// Mean and standard error of max(payoff, 0) over antithetic pairs.
CheyettePrice optionValue(const State& s, const Payoff& payoff) {
    const std::size_t half = s.size() / 2;
    double sum = 0.0, sumSq = 0.0;
    for (std::size_t p = 0; p < half; ++p) {
        const double a = std::max(payoff(s.x[p], s.y[p], s.integral[p]), 0.0), b = std::max(payoff(s.x[p + half], s.y[p + half], s.integral[p + half]), 0.0);
        const double v = 0.5 * (a + b);
        sum += v;
        sumSq += v * v;
    }
    const double n = static_cast<double>(half), mean = sum / n;
    return {mean, std::sqrt(std::max(sumSq / n - mean * mean, 0.0) / n)};
}

void check(const CheyetteOptions& o, const char* who) {
    RL_REQUIRE(o.paths >= 200 && o.paths % 2 == 0, std::string(who) + ": need an even number of at least 200 paths");
    RL_REQUIRE(o.stepsPerYear >= 1, std::string(who) + ": need at least one time step a year");
}

/// Walk a state through the model's own level intervals from t0 to t1, drawing normals as it goes.
void advanceThrough(State& s, const Cheyette& model, double t0, double t1, int stepsPerYear, math::Rng& rng) {
    double from = t0;
    while (from < t1 - 1e-12) {
        double to = t1;
        for (double node : model.times())
            if (node > from + 1e-12 && node < to - 1e-12) to = node;
        const int steps = stepsBetween(from, to, stepsPerYear);
        advance(s, model, model.level(to), from, to, steps, drawNormals(rng, s.size() / 2, steps));
        from = to;
    }
}

} // namespace

std::vector<CheyettePrice> priceEuropeansCheyette(const Cheyette& model, const std::vector<models::ExerciseRight>& rights, const CheyetteOptions& opt) {
    check(opt, "priceEuropeansCheyette");
    for (std::size_t i = 1; i < rights.size(); ++i) RL_REQUIRE(rights[i].time >= rights[i - 1].time, "priceEuropeansCheyette: rights must be sorted by exercise time");
    State s(opt.paths);
    math::Rng rng(opt.seed);
    std::vector<CheyettePrice> out;
    double t = 0.0;
    for (const models::ExerciseRight& r : rights) {
        RL_REQUIRE(r.time > 0.0, "priceEuropeansCheyette: exercise times must be positive");
        if (r.time > t + 1e-12) {
            advanceThrough(s, model, t, r.time, opt.stepsPerYear, rng);
            t = r.time;
        }
        out.push_back(optionValue(s, Payoff(model, r)));
    }
    return out;
}

CheyetteCalibration calibrateCheyette(double kappa, double skew, const std::vector<models::CalibrationInstrument>& strip, const CheyetteOptions& opt) {
    check(opt, "calibrateCheyette");
    RL_REQUIRE(!strip.empty(), "calibrateCheyette: empty strip");
    for (std::size_t i = 1; i < strip.size(); ++i)
        RL_REQUIRE(strip[i].right.time > strip[i - 1].right.time, "calibrateCheyette: the strip must be in increasing expiry order");

    State state(opt.paths);
    math::Rng rng(opt.seed);
    std::vector<double> times, levels, prices;
    double t = 0.0;
    for (const models::CalibrationInstrument& ci : strip) {
        const double end = ci.right.time;
        const int steps = stepsBetween(t, end, opt.stepsPerYear);
        const std::vector<double> normals = drawNormals(rng, opt.paths / 2, steps);
        times.push_back(end);
        levels.push_back(0.01);
        // G does not depend on the levels, so one model object serves every trial.
        const Cheyette shape(kappa, skew, times, levels);
        const Payoff payoff(shape, ci.right);
        State trial(0);
        auto error = [&](double logLevel) {
            trial = state;
            advance(trial, shape, std::exp(logLevel), t, end, steps, normals);
            return optionValue(trial, payoff).value - ci.marketPrice;
        };
        double lo = std::log(1e-5), hi = std::log(0.004);
        if (error(lo) >= 0.0)
            throw CalibrationError("calibrateCheyette: " + ci.right.label + " needs less variance than the earlier expiries already imply", {error(lo)});
        while (error(hi) <= 0.0) {
            hi += std::log(2.0);
            if (hi > std::log(0.5)) throw CalibrationError("calibrateCheyette: market price of " + ci.right.label + " is not attainable", {error(hi)});
        }
        math::RootOptions ro;
        ro.xTol = 1e-11;
        ro.maxIter = 100;
        levels.back() = std::exp(math::brent(error, lo, hi, ro).x);
        advance(state, shape, levels.back(), t, end, steps, normals);
        prices.push_back(optionValue(state, payoff).value);
        t = end;
    }

    CheyetteCalibration out{Cheyette(kappa, skew, times, levels), {}, 0.0, true};
    for (std::size_t i = 0; i < strip.size(); ++i) {
        const models::ExerciseRight& r = strip[i].right;
        CheyetteCalibrationRow row{r.label, r.time, strip[i].marketVol * 1e4, models::impliedNormalVol(r, prices[i]) * 1e4, 0.0, levels[i] * 1e4};
        row.errorBp = row.modelVolBp - row.marketVolBp;
        out.maxAbsErrorBp = std::max(out.maxAbsErrorBp, std::fabs(row.errorBp));
        out.rows.push_back(std::move(row));
    }
    out.ok = out.maxAbsErrorBp < 0.01;
    return out;
}

CheyetteLsmcResult priceBermudanCheyette(const Cheyette& model, const std::vector<models::ExerciseRight>& rights, const CheyetteLsmcOptions& opt) {
    const std::size_t n = rights.size();
    RL_REQUIRE(n > 0, "priceBermudanCheyette: no exercise dates");
    for (std::size_t i = 1; i < n; ++i) RL_REQUIRE(rights[i].time > rights[i - 1].time, "priceBermudanCheyette: rights must be in increasing exercise time");
    check(opt.pricing, "priceBermudanCheyette");
    RL_REQUIRE(opt.trainingPaths >= 1000 && opt.trainingPaths % 2 == 0, "priceBermudanCheyette: need an even number of at least 1000 training paths");
    RL_REQUIRE(opt.degree >= 1 && opt.degree <= 6, "priceBermudanCheyette: degree must be in [1, 6]");
    std::vector<Payoff> payoffs;
    for (const auto& r : rights) payoffs.emplace_back(model, r);

    // Regressors: powers of the standardised short-rate state, plus the standardised accumulated
    // variance and its product with the state. With no skew y is the same on every path and drops out.
    struct Scale {
        double mx = 0.0, sx = 1.0, my = 0.0, sy = 0.0;
    };
    std::vector<Scale> scales(n);
    auto basis = [&](std::size_t i, double x, double y, std::vector<double>& row) {
        const Scale& s = scales[i];
        const double u = (x - s.mx) / s.sx;
        row.clear();
        double pu = 1.0;
        for (int a = 0; a <= opt.degree; ++a, pu *= u) row.push_back(pu);
        if (s.sy > 0.0) {
            const double v = (y - s.my) / s.sy;
            row.push_back(v);
            row.push_back(u * v);
        }
    };

    // --- Pass 1: fit the rule -------------------------------------------------------------------------
    const std::size_t N = opt.trainingPaths;
    std::vector<std::vector<double>> xs(n), ys(n), hs(n), discount(n);
    {
        State s(N);
        math::Rng rng(math::SeedPolicy::derive(opt.pricing.seed, 7));
        double t = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            advanceThrough(s, model, t, rights[i].time, opt.pricing.stepsPerYear, rng);
            t = rights[i].time;
            xs[i] = s.x;
            ys[i] = s.y;
            discount[i].resize(N);
            for (std::size_t p = 0; p < N; ++p) discount[i][p] = std::exp(-s.integral[p]);
            hs[i].resize(N);
            double mx = 0.0, my = 0.0, vx = 0.0, vy = 0.0;
            for (std::size_t p = 0; p < N; ++p) {
                hs[i][p] = payoffs[i](s.x[p], s.y[p], s.integral[p]);
                mx += s.x[p];
                my += s.y[p];
            }
            mx /= static_cast<double>(N);
            my /= static_cast<double>(N);
            for (std::size_t p = 0; p < N; ++p) {
                vx += (s.x[p] - mx) * (s.x[p] - mx);
                vy += (s.y[p] - my) * (s.y[p] - my);
            }
            const double sx = std::sqrt(vx / static_cast<double>(N)), sy = std::sqrt(vy / static_cast<double>(N));
            scales[i] = {mx, sx > 0.0 ? sx : 1.0, my, sy > 1e-9 * std::fabs(my) ? sy : 0.0};
        }
    }
    std::vector<std::vector<double>> coefficients(n);
    std::vector<double> row;
    auto continuation = [&](std::size_t i, double x, double y) {
        basis(i, x, y, row);
        double c = 0.0;
        for (std::size_t k = 0; k < row.size(); ++k) c += coefficients[i][k] * row[k];
        return c;
    };
    // Exercise value and continuation are both compared as values AT the decision date.
    auto exercise = [&](std::size_t i, double x, double y, double h) {
        if (h <= 0.0) return false;
        if (i + 1 == n) return true;
        if (coefficients[i].empty()) return false;
        return h >= continuation(i, x, y);
    };
    // Cash flows are discounted to today. The regression explains them by the state at date i, so
    // they are first expressed in units of the discount factor to date i (i.e. as values at date i).
    std::vector<double> cash(N);
    for (std::size_t p = 0; p < N; ++p) cash[p] = std::max(hs[n - 1][p], 0.0);
    for (std::size_t i = n - 1; i-- > 0;) {
        std::vector<std::size_t> itm;
        for (std::size_t p = 0; p < N; ++p)
            if (hs[i][p] > 0.0) itm.push_back(p);
        basis(i, 0.0, 0.0, row);
        const std::size_t terms = row.size();
        if (itm.size() < 10 * terms) continue;
        Eigen::MatrixXd A(static_cast<Eigen::Index>(itm.size()), static_cast<Eigen::Index>(terms));
        Eigen::VectorXd b(static_cast<Eigen::Index>(itm.size()));
        for (std::size_t r = 0; r < itm.size(); ++r) {
            const std::size_t p = itm[r];
            basis(i, xs[i][p], ys[i][p], row);
            for (std::size_t k = 0; k < terms; ++k) A(static_cast<Eigen::Index>(r), static_cast<Eigen::Index>(k)) = row[k];
            b[static_cast<Eigen::Index>(r)] = cash[p] / discount[i][p];
        }
        const Eigen::VectorXd fit = A.colPivHouseholderQr().solve(b);
        coefficients[i].assign(fit.data(), fit.data() + fit.size());
        for (std::size_t p : itm)
            if (hs[i][p] / discount[i][p] >= continuation(i, xs[i][p], ys[i][p])) cash[p] = hs[i][p];
    }

    // --- Pass 2: value the rule on fresh antithetic paths ------------------------------------------------
    const std::size_t M = opt.pricing.paths, half = M / 2;
    State s(M);
    math::Rng rng(opt.pricing.seed);
    std::vector<double> value(M, 0.0);
    std::vector<char> done(M, 0);
    std::vector<std::size_t> exercised(n, 0);
    double t = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        advanceThrough(s, model, t, rights[i].time, opt.pricing.stepsPerYear, rng);
        t = rights[i].time;
        for (std::size_t p = 0; p < M; ++p) {
            if (done[p]) continue;
            const double h = payoffs[i](s.x[p], s.y[p], s.integral[p]);
            if (exercise(i, s.x[p], s.y[p], h * std::exp(s.integral[p]))) {
                value[p] = h;
                done[p] = 1;
                ++exercised[i];
            }
        }
    }
    double sum = 0.0, sumSq = 0.0;
    for (std::size_t p = 0; p < half; ++p) {
        const double v = 0.5 * (value[p] + value[p + half]);
        sum += v;
        sumSq += v * v;
    }
    CheyetteLsmcResult res;
    const double pairs = static_cast<double>(half);
    res.lower = sum / pairs;
    res.standardError = std::sqrt(std::max(sumSq / pairs - res.lower * res.lower, 0.0) / pairs);
    for (std::size_t c : exercised) res.exerciseProbability.push_back(static_cast<double>(c) / static_cast<double>(M));
    return res;
}

} // namespace rl::pricers
