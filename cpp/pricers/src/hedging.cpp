#include "rl/pricers/hedging.hpp"

#include "rl/core/error.hpp"
#include "rl/math/rng.hpp"

#include <algorithm>
#include <cmath>

namespace rl::pricers {

namespace {

const double kTwoPi = 6.283185307179586;

/// The Bermudan's value on the first exercise date, on a grid, by backward induction with exact
/// Gaussian transitions (the method of xva/bermudan_exposure.cpp).
struct FirstDateValue {
    double lo = 0.0, step = 1.0;
    std::vector<double> v; ///< max(exercise, continuation) at t_1
};

FirstDateValue firstDateValue(const BermudanProblem& p, int m) {
    const std::size_t n = p.size(), points = static_cast<std::size_t>(m);
    std::vector<double> next, grid(points);
    double nextLo = 0.0, nextStep = 1.0;
    FirstDateValue out;
    for (std::size_t i = n; i-- > 0;) {
        const double sd = std::sqrt(p.zeta(i)), lo = -8.0 * sd, step = 16.0 * sd / static_cast<double>(points - 1);
        std::vector<double> value(points);
        for (std::size_t a = 0; a < points; ++a) {
            const double x = lo + step * static_cast<double>(a);
            double cont = 0.0;
            if (i + 1 < n) {
                const double var = p.zeta(i + 1) - p.zeta(i), norm = nextStep / std::sqrt(kTwoPi * var);
                for (std::size_t b = 0; b < points; ++b) {
                    const double d = nextLo + nextStep * static_cast<double>(b) - x;
                    cont += next[b] * std::exp(-0.5 * d * d / var);
                }
                cont *= norm;
            }
            value[a] = std::max(p.exerciseValue(i, x), cont);
        }
        next.swap(value);
        nextLo = lo;
        nextStep = step;
    }
    out.lo = nextLo;
    out.step = nextStep;
    out.v = std::move(next);
    return out;
}

/// V(t, x) and ∂V/∂x for t before the first exercise date, from the value there and the remaining variance.
std::pair<double, double> valueAndDelta(const FirstDateValue& f, double remainingVariance, double x) {
    if (remainingVariance <= 1e-14) {
        // On the date itself: linear interpolation on the grid, and its slope.
        const double u = std::clamp((x - f.lo) / f.step, 0.0, static_cast<double>(f.v.size() - 1) - 1e-9);
        const std::size_t i = static_cast<std::size_t>(u);
        const double w = u - static_cast<double>(i);
        return {(1.0 - w) * f.v[i] + w * f.v[i + 1], (f.v[i + 1] - f.v[i]) / f.step};
    }
    const double norm = f.step / std::sqrt(kTwoPi * remainingVariance);
    double v = 0.0, dv = 0.0;
    for (std::size_t b = 0; b < f.v.size(); ++b) {
        const double d = f.lo + f.step * static_cast<double>(b) - x, k = std::exp(-0.5 * d * d / remainingVariance);
        v += f.v[b] * k;
        dv += f.v[b] * k * d / remainingVariance;
    }
    return {v * norm, dv * norm};
}

/// The hedge swap: the swap delivered on the first exercise date, deflated: Σ c_k exp(−H_k x − ½H_k² ζ(t)).
struct HedgeSwap {
    std::vector<double> h, c;
    HedgeSwap(const models::Lgm1F& model, const models::ExerciseRight& r) {
        for (const auto& f : r.flows) {
            h.push_back(model.H(f.time));
            c.push_back(f.amount * f.df);
        }
    }
    std::pair<double, double> operator()(double zeta, double x) const noexcept {
        double v = 0.0, dv = 0.0;
        for (std::size_t k = 0; k < c.size(); ++k) {
            const double term = c[k] * std::exp(-h[k] * x - 0.5 * h[k] * h[k] * zeta);
            v += term;
            dv -= h[k] * term;
        }
        return {v, dv};
    }
};

models::Lgm1F scaled(const models::Lgm1F& m, double s) {
    std::vector<double> z = m.zetaValues();
    for (double& v : z) v *= s * s;
    return models::Lgm1F(m.kappa(), m.zetaTimes(), z);
}

} // namespace

HedgingResult simulateDeltaHedge(const models::Lgm1F& trueModel, const std::vector<models::ExerciseRight>& rights, const HedgingOptions& opt) {
    RL_REQUIRE(!rights.empty(), "simulateDeltaHedge: no exercise rights");
    RL_REQUIRE(opt.paths >= 100 && opt.rebalancesPerYear >= 1 && opt.volScale > 0.0, "simulateDeltaHedge: need at least 100 paths, one rebalance a year and a positive vol scale");
    RL_REQUIRE(opt.gridPoints >= 101 && opt.gridPoints % 2 == 1, "simulateDeltaHedge: need an odd grid of at least 101 points");
    const models::Lgm1F hedgerModel = scaled(trueModel, opt.volScale);
    const BermudanProblem truth(trueModel, rights), hedger(hedgerModel, rights);
    const FirstDateValue fTrue = firstDateValue(truth, opt.gridPoints), fHedger = firstDateValue(hedger, opt.gridPoints);
    const HedgeSwap swapTrue(trueModel, rights.front()), swapHedger(hedgerModel, rights.front());
    const double T = rights.front().time, zT = trueModel.zeta(T), zH = hedgerModel.zeta(T);

    HedgingResult out;
    out.horizon = T;
    out.truePrice = valueAndDelta(fTrue, zT, 0.0).first;
    out.premium = valueAndDelta(fHedger, zH, 0.0).first;

    const int steps = std::max(1, static_cast<int>(std::lround(T * opt.rebalancesPerYear)));
    std::vector<double> times(static_cast<std::size_t>(steps) + 1);
    for (int j = 0; j <= steps; ++j) times[static_cast<std::size_t>(j)] = T * j / steps;
    // The hedger's delta at each rebalance depends only on (t, x): precompute nothing, the integral is cheap enough.
    math::Rng rng(opt.seed);
    double sum = 0.0, sumSq = 0.0;
    for (std::size_t p = 0; p < opt.paths; ++p) {
        double x = 0.0, pnl = out.premium;
        for (int j = 0; j < steps; ++j) {
            const double t0 = times[static_cast<std::size_t>(j)], t1 = times[static_cast<std::size_t>(j) + 1];
            const double z0 = trueModel.zeta(t0), z1 = trueModel.zeta(t1);
            double hedgeRatio = 0.0;
            if (opt.hedge) {
                const double dv = valueAndDelta(fHedger, zH - hedgerModel.zeta(t0), x).second;
                const double ds = swapHedger(hedgerModel.zeta(t0), x).second;
                hedgeRatio = dv / ds;
            }
            const double h0 = swapTrue(z0, x).first;
            x += std::sqrt(std::max(z1 - z0, 0.0)) * rng.gaussian();
            // Short the option, long hedgeRatio swaps: the gain is the swap's change.
            pnl += hedgeRatio * (swapTrue(z1, x).first - h0);
        }
        pnl -= valueAndDelta(fTrue, 0.0, x).first; // buy back at the market's value
        sum += pnl;
        sumSq += pnl * pnl;
    }
    const double n = static_cast<double>(opt.paths);
    out.meanPnl = sum / n;
    out.pnlStdDev = std::sqrt(std::max(sumSq / n - out.meanPnl * out.meanPnl, 0.0));
    out.meanPnlStandardError = out.pnlStdDev / std::sqrt(n);
    return out;
}

namespace {
double onebp(const bootstrap::BootstrapInstrument& bi) { return bi.type == market::QuoteType::Future ? -0.01 : 1e-4; }
} // namespace

models::UsdMarket applyMove(const models::UsdMarket& market, const MarketMove& move) {
    const std::size_t n = market.problem.size();
    RL_REQUIRE(move.quoteShiftBp.empty() || move.quoteShiftBp.size() == n, "applyMove: need one shift per curve quote");
    math::Vector q = market.problem.quoteValues();
    bool curve = false;
    for (std::size_t i = 0; i < move.quoteShiftBp.size(); ++i)
        if (move.quoteShiftBp[i] != 0.0) {
            q[static_cast<Eigen::Index>(i)] += move.quoteShiftBp[i] * onebp(market.problem.instruments()[i]);
            curve = true;
        }
    models::UsdMarket out = curve ? market.withCurveQuotes(q) : market;
    if (move.volShiftBp != 0.0) out = out.withVolShift(move.volShiftBp * 1e-4);
    return out;
}

Sensitivities computeSensitivities(const models::UsdMarket& market, const std::function<double(const models::UsdMarket&)>& price, double gammaBumpBp) {
    RL_REQUIRE(gammaBumpBp > 0.0, "computeSensitivities: the gamma bump must be positive");
    const std::size_t n = market.problem.size();
    auto shifted = [&](std::size_t i, double bp) {
        MarketMove m;
        m.quoteShiftBp.assign(n, 0.0);
        if (i < n) m.quoteShiftBp[i] = bp;
        else m.quoteShiftBp.assign(n, bp);
        return price(applyMove(market, m));
    };
    Sensitivities s;
    s.value = price(market);
    s.delta.resize(n);
    for (std::size_t i = 0; i < n; ++i) s.delta[i] = 0.5 * (shifted(i, 1.0) - shifted(i, -1.0));
    s.parallelGamma = (shifted(n, gammaBumpBp) - 2.0 * s.value + shifted(n, -gammaBumpBp)) / (gammaBumpBp * gammaBumpBp);
    MarketMove vu, vd;
    vu.volShiftBp = 1.0;
    vd.volShiftBp = -1.0;
    s.vega = 0.5 * (price(applyMove(market, vu)) - price(applyMove(market, vd)));
    return s;
}

PnlExplain explainPnl(const models::UsdMarket& market, const Sensitivities& s, const std::function<double(const models::UsdMarket&)>& price, const MarketMove& move) {
    PnlExplain e;
    e.actual = price(applyMove(market, move)) - s.value;
    double mean = 0.0;
    for (std::size_t i = 0; i < move.quoteShiftBp.size(); ++i) {
        e.delta += s.delta[i] * move.quoteShiftBp[i];
        mean += move.quoteShiftBp[i];
    }
    if (!move.quoteShiftBp.empty()) mean /= static_cast<double>(move.quoteShiftBp.size());
    e.gamma = 0.5 * s.parallelGamma * mean * mean;
    e.vega = s.vega * move.volShiftBp;
    return e;
}

} // namespace rl::pricers
