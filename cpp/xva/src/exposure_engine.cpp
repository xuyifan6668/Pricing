#include "rl/xva/exposure_engine.hpp"

#include "rl/core/error.hpp"
#include "rl/math/rng.hpp"
#include "rl/models/g2pp.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <variant>

namespace rl::xva {

namespace {

using instruments::FixedCoupon;
using instruments::OvernightCoupon;

const double kEps = 1e-9; // time comparisons: a millisecond of a year

// The engine is written once for a two-dimensional driftless Gaussian state. A deflated bond maturing
// at T is D(T) exp(−L(T)·x − ½ L(T)ᵀ C(t) L(T)); the one-factor model is the case where the second
// loading and the second state are identically zero.
struct Load {
    double h1 = 0.0, h2 = 0.0;
};
struct Cov {
    double c11 = 0.0, c12 = 0.0, c22 = 0.0;
};
double quad(const Load& l, const Cov& c) noexcept { return l.h1 * l.h1 * c.c11 + 2.0 * l.h1 * l.h2 * c.c12 + l.h2 * l.h2 * c.c22; }

struct Dynamics {
    bool twoFactor = false;
    std::function<Load(double)> loading;   ///< L(T)
    std::function<Cov(double)> covariance; ///< C(t)
};

struct FixedItem {
    double end = 0.0;   ///< accrual end: the coupon is in the valuation until then
    Load h;             ///< L(payment)
    double c = 0.0;     ///< sign · N · K · τ · D(payment)
};
struct FloatItem {
    double start = 0.0, end = 0.0;
    Load hs, he;
    double ds = 0.0, de = 0.0; ///< D(start), D(end)
    double n = 0.0;            ///< sign · N · D(payment)/D(end): notional with the payment delay folded in
    std::size_t id = 0;        ///< index into the per-path array of fixed rates
};
struct Trade {
    std::vector<FixedItem> fixed;
    std::vector<FloatItem> floating;
};

Trade compile(const Dynamics& model, const instruments::Swap& swap, const market::CurveSet<double>& curves, std::size_t& nextId) {
    const market::YieldCurve<double>& disc = curves.discount();
    Trade t;
    for (const auto& cf : swap.fixedLeg().cashflows()) {
        const auto* c = std::get_if<FixedCoupon>(&cf);
        RL_REQUIRE(c != nullptr, "simulateExposure: the fixed leg of " + swap.id() + " must hold fixed coupons only");
        const double pay = disc.time(c->payment);
        t.fixed.push_back({disc.time(c->accrualEnd), model.loading(pay), swap.fixedLeg().sign() * c->notional * c->rate * c->tau * disc.df(pay)});
    }
    for (const auto& cf : swap.floatLeg().cashflows()) {
        const auto* c = std::get_if<OvernightCoupon>(&cf);
        RL_REQUIRE(c != nullptr, "simulateExposure: the floating leg of " + swap.id() + " must hold overnight coupons only");
        RL_REQUIRE(c->spread == 0.0, "simulateExposure: a spread on the floating leg is not supported");
        RL_REQUIRE(&curves.projection(c->index.name) == &disc, "simulateExposure: " + c->index.name + " is not projected off the discount curve (single-curve model)");
        FloatItem f;
        f.start = disc.time(c->accrualStart);
        f.end = disc.time(c->accrualEnd);
        f.hs = model.loading(std::max(f.start, 0.0));
        f.he = model.loading(f.end);
        f.ds = disc.df(std::max(f.start, 0.0));
        f.de = disc.df(f.end);
        f.n = swap.floatLeg().sign() * c->notional * disc.df(c->payment) / f.de;
        f.id = nextId++;
        t.floating.push_back(f);
    }
    return t;
}

/// What a trade is worth at one valuation date, as a function of the state: bonds still to come,
/// plus the floating coupons whose rate was fixed earlier on the path.
struct Term {
    double h1, h2, c;
};
struct InProgress {
    double h1, h2, c;
    std::size_t id;
};
struct Valuation {
    std::vector<std::vector<Term>> terms;        // per trade
    std::vector<std::vector<InProgress>> fixedRates; // per trade
};

Valuation valuationAt(const std::vector<Trade>& trades, double t, const Cov& cov) {
    Valuation v;
    for (const Trade& tr : trades) {
        std::vector<Term> terms;
        std::vector<InProgress> inProgress;
        auto bond = [&](const Load& h, double c) { terms.push_back({h.h1, h.h2, c * std::exp(-0.5 * quad(h, cov))}); };
        for (const FixedItem& f : tr.fixed)
            if (f.end > t + kEps) bond(f.h, f.c);
        for (const FloatItem& f : tr.floating) {
            if (f.end <= t + kEps) continue;
            if (f.start >= t - kEps) {
                bond(f.hs, f.n * f.ds);
                bond(f.he, -f.n * f.de);
            } else {
                inProgress.push_back({f.he.h1, f.he.h2, f.n * f.de * std::exp(-0.5 * quad(f.he, cov)), f.id});
            }
        }
        v.terms.push_back(std::move(terms));
        v.fixedRates.push_back(std::move(inProgress));
    }
    return v;
}

double tradeValue(const Valuation& v, std::size_t trade, double x1, double x2, const std::vector<double>& rates) {
    double value = 0.0;
    for (const Term& t : v.terms[trade]) value += t.c * std::exp(-t.h1 * x1 - t.h2 * x2);
    for (const InProgress& p : v.fixedRates[trade]) value += p.c * std::exp(-p.h1 * x1 - p.h2 * x2) * rates[p.id];
    return value;
}

/// The collateral the agreement calls for against a portfolio worth `value`.
double collateralCall(const Csa& csa, double value) {
    double call = 0.0;
    if (value > csa.threshold) call = value - csa.threshold;
    else if (value < -csa.threshold) call = value + csa.threshold;
    return std::fabs(call) >= csa.minimumTransfer ? call : 0.0;
}

} // namespace

std::vector<double> periodStartTimes(const instruments::Swap& swap, const market::CurveSet<double>& curves) {
    std::vector<double> out{0.0};
    for (const auto& cf : swap.fixedLeg().cashflows()) {
        const double t = curves.discount().time(std::get<FixedCoupon>(cf).accrualStart);
        if (t > out.back() + kEps) out.push_back(t);
    }
    return out;
}

std::vector<double> exposureTimes(const std::vector<instruments::Swap>& trades, const market::CurveSet<double>& curves, double stepYears) {
    RL_REQUIRE(stepYears > 0.0, "exposureTimes: the step must be positive");
    std::vector<double> all{0.0};
    double horizon = 0.0;
    for (const auto& s : trades) {
        for (const auto& cf : s.fixedLeg().cashflows()) {
            const auto& c = std::get<FixedCoupon>(cf);
            all.push_back(curves.discount().time(c.accrualStart));
            horizon = std::max(horizon, curves.discount().time(c.accrualEnd));
        }
    }
    for (double t = stepYears; t < horizon - kEps; t += stepYears) all.push_back(t);
    std::sort(all.begin(), all.end());
    std::vector<double> out;
    for (double t : all)
        if (t >= 0.0 && t < horizon - kEps && (out.empty() || t > out.back() + 1e-6)) out.push_back(t);
    return out;
}

namespace {

ExposureProfile simulate(const Dynamics& model, const std::vector<instruments::Swap>& swaps, const market::CurveSet<double>& curves,
                         const std::vector<double>& times, const std::optional<Csa>& csa, const ExposureOptions& opt) {
    RL_REQUIRE(!swaps.empty(), "simulateExposure: no trades");
    RL_REQUIRE(!times.empty() && times.front() == 0.0, "simulateExposure: exposure times must start at 0");
    for (std::size_t i = 1; i < times.size(); ++i) RL_REQUIRE(times[i] > times[i - 1] + kEps, "simulateExposure: exposure times must be strictly increasing");
    RL_REQUIRE(opt.paths >= 100 && opt.quantile > 0.5 && opt.quantile < 1.0, "simulateExposure: need at least 100 paths and a quantile in (0.5, 1)");
    if (csa) RL_REQUIRE(csa->threshold >= 0.0 && csa->minimumTransfer >= 0.0 && csa->marginPeriodYears >= 0.0, "simulateExposure: CSA terms must be non-negative");
    const market::YieldCurve<double>& disc = curves.discount();

    std::size_t rateCount = 0;
    std::vector<Trade> trades;
    for (const auto& s : swaps) trades.push_back(compile(model, s, curves, rateCount));
    const std::size_t nTrades = trades.size(), nTimes = times.size();

    // --- The simulation dates: exposure dates, the dates collateral is called on, and the dates on
    //     which a floating coupon's rate is fixed. --------------------------------------------------
    struct Reset {
        Load hs, he;
        double ratio; ///< D(start)/D(end)
        std::size_t id;
    };
    struct Step {
        double time = 0.0;
        Cov cov;                         ///< C(t)
        double a11 = 0.0, a21 = 0.0, a22 = 0.0; ///< Cholesky factor of C(t) − C(previous step)
        Load h;                          ///< L(t): for the numeraire
        double df = 1.0;                 ///< D(t)
        std::ptrdiff_t exposure = -1;   ///< index into `times`, or −1
        std::ptrdiff_t collateralFor = -1; ///< the exposure date whose collateral is called here, or −1
        bool valued = false;
        Valuation valuation;
        std::vector<Reset> resets;
    };
    std::vector<double> stepTimes;
    for (std::size_t i = 1; i < nTimes; ++i) stepTimes.push_back(times[i]);
    if (csa)
        for (std::size_t i = 1; i < nTimes; ++i)
            if (times[i] - csa->marginPeriodYears > kEps) stepTimes.push_back(times[i] - csa->marginPeriodYears);
    for (const Trade& t : trades)
        for (const FloatItem& f : t.floating)
            if (f.start > kEps && f.start < times.back()) stepTimes.push_back(f.start);
    std::sort(stepTimes.begin(), stepTimes.end());
    std::vector<Step> steps;
    for (double t : stepTimes) {
        if (!steps.empty() && t < steps.back().time + kEps) continue;
        Step s;
        s.time = t;
        s.cov = model.covariance(t);
        const Cov before = steps.empty() ? Cov{} : steps.back().cov;
        const double d11 = std::max(s.cov.c11 - before.c11, 0.0), d12 = s.cov.c12 - before.c12, d22 = std::max(s.cov.c22 - before.c22, 0.0);
        s.a11 = std::sqrt(d11);
        s.a21 = s.a11 > 0.0 ? d12 / s.a11 : 0.0;
        s.a22 = std::sqrt(std::max(d22 - s.a21 * s.a21, 0.0));
        s.h = model.loading(t);
        s.df = disc.df(t);
        steps.push_back(std::move(s));
    }
    auto stepAt = [&](double t) -> Step& {
        for (Step& s : steps)
            if (std::fabs(s.time - t) < kEps) return s;
        throw NumericalError("simulateExposure: a simulation date went missing");
    };
    for (std::size_t i = 1; i < nTimes; ++i) {
        Step& s = stepAt(times[i]);
        s.exposure = static_cast<std::ptrdiff_t>(i);
        s.valued = true;
        if (csa && times[i] - csa->marginPeriodYears > kEps) {
            Step& c = stepAt(times[i] - csa->marginPeriodYears);
            c.collateralFor = static_cast<std::ptrdiff_t>(i);
            c.valued = true;
        }
    }
    for (const Trade& t : trades)
        for (const FloatItem& f : t.floating)
            if (f.start > kEps && f.start < times.back()) stepAt(f.start).resets.push_back({f.hs, f.he, f.ds / f.de, f.id});
    for (Step& s : steps)
        if (s.valued) s.valuation = valuationAt(trades, s.time, s.cov);

    // Coupons already accruing today have their rate fixed off today's curve.
    std::vector<double> ratesToday(rateCount, 0.0);
    for (const Trade& t : trades)
        for (const FloatItem& f : t.floating)
            if (f.start <= kEps) ratesToday[f.id] = f.ds / f.de - 1.0;

    ExposureProfile out;
    out.times = times;
    out.epe.assign(nTimes, 0.0);
    out.ene.assign(nTimes, 0.0);
    out.epeStandardError.assign(nTimes, 0.0);
    out.pfe.assign(nTimes, 0.0);
    out.grossEpe.assign(nTimes, 0.0);
    if (csa) {
        out.collateralisedEpe.assign(nTimes, 0.0);
        out.collateralisedEne.assign(nTimes, 0.0);
    }

    // --- Today ---------------------------------------------------------------------------------------
    const Valuation today = valuationAt(trades, 0.0, Cov{});
    double collateralToday = 0.0;
    for (std::size_t k = 0; k < nTrades; ++k) {
        const double v = tradeValue(today, k, 0.0, 0.0, ratesToday);
        out.valueToday += v;
        out.grossEpe[0] += std::max(v, 0.0);
    }
    out.epe[0] = out.pfe[0] = std::max(out.valueToday, 0.0);
    out.ene[0] = std::min(out.valueToday, 0.0);
    if (csa) {
        collateralToday = collateralCall(*csa, out.valueToday);
        out.collateralisedEpe[0] = std::max(out.valueToday - collateralToday, 0.0);
        out.collateralisedEne[0] = std::min(out.valueToday - collateralToday, 0.0);
    }

    // --- Paths ---------------------------------------------------------------------------------------
    std::vector<double> sumSq(nTimes, 0.0);
    std::vector<std::vector<double>> positive(nTimes, std::vector<double>(opt.paths, 0.0)); // for the quantile
    std::vector<double> rates(rateCount), collateral(nTimes);
    math::Rng rng(opt.seed);
    for (std::size_t p = 0; p < opt.paths; ++p) {
        rates = ratesToday;
        std::fill(collateral.begin(), collateral.end(), collateralToday); // margin period reaching back to today
        double x1 = 0.0, x2 = 0.0;
        for (const Step& s : steps) {
            const double z1 = rng.gaussian();
            x1 += s.a11 * z1;
            if (model.twoFactor) x2 += s.a21 * z1 + s.a22 * rng.gaussian();
            // P(s, e; x) = D(e)/D(s) exp(−(L_e − L_s)·x − ½ (L_eᵀ C L_e − L_sᵀ C L_s)): the term rate fixed here.
            for (const Reset& r : s.resets)
                rates[r.id] = r.ratio * std::exp((r.he.h1 - r.hs.h1) * x1 + (r.he.h2 - r.hs.h2) * x2 + 0.5 * (quad(r.he, s.cov) - quad(r.hs, s.cov))) - 1.0;
            if (!s.valued) continue;
            double net = 0.0, gross = 0.0;
            for (std::size_t k = 0; k < nTrades; ++k) {
                const double v = tradeValue(s.valuation, k, x1, x2, rates);
                net += v;
                gross += std::max(v, 0.0);
            }
            if (s.collateralFor >= 0) {
                // The agreement works in money: undo the deflation, apply the terms, deflate again.
                const double numeraire = std::exp(s.h.h1 * x1 + s.h.h2 * x2 + 0.5 * quad(s.h, s.cov)) / s.df;
                collateral[static_cast<std::size_t>(s.collateralFor)] = collateralCall(*csa, net * numeraire) / numeraire;
            }
            if (s.exposure < 0) continue;
            const std::size_t i = static_cast<std::size_t>(s.exposure);
            const double pos = std::max(net, 0.0);
            out.epe[i] += pos;
            sumSq[i] += pos * pos;
            out.ene[i] += std::min(net, 0.0);
            out.grossEpe[i] += gross;
            positive[i][p] = pos;
            if (csa) {
                const double residual = net - collateral[i];
                out.collateralisedEpe[i] += std::max(residual, 0.0);
                out.collateralisedEne[i] += std::min(residual, 0.0);
            }
        }
    }
    const double n = static_cast<double>(opt.paths);
    const std::size_t rank = std::min(opt.paths - 1, static_cast<std::size_t>(opt.quantile * n));
    for (std::size_t i = 1; i < nTimes; ++i) {
        out.epe[i] /= n;
        out.ene[i] /= n;
        out.grossEpe[i] /= n;
        out.epeStandardError[i] = std::sqrt(std::max(sumSq[i] / n - out.epe[i] * out.epe[i], 0.0) / n);
        std::nth_element(positive[i].begin(), positive[i].begin() + static_cast<std::ptrdiff_t>(rank), positive[i].end());
        out.pfe[i] = positive[i][rank];
        if (csa) {
            out.collateralisedEpe[i] /= n;
            out.collateralisedEne[i] /= n;
        }
    }
    return out;
}

/// The latest date any trade's cashflows reach: the two-factor model's loadings are measured from it.
double horizonOf(const std::vector<instruments::Swap>& swaps, const market::CurveSet<double>& curves) {
    double h = 0.0;
    for (const auto& s : swaps)
        for (const auto& cf : s.fixedLeg().cashflows()) h = std::max(h, curves.discount().time(std::get<FixedCoupon>(cf).payment));
    return h;
}

} // namespace

ExposureProfile simulateExposure(const models::Lgm1F& model, const std::vector<instruments::Swap>& swaps, const market::CurveSet<double>& curves,
                                 const std::vector<double>& times, const std::optional<Csa>& csa, const ExposureOptions& opt) {
    Dynamics d;
    d.twoFactor = false;
    d.loading = [&](double T) { return Load{model.H(T), 0.0}; };
    d.covariance = [&](double t) { return Cov{model.zeta(t), 0.0, 0.0}; };
    return simulate(d, swaps, curves, times, csa, opt);
}

ExposureProfile simulateExposure(const models::G2pp& model, const std::vector<instruments::Swap>& swaps, const market::CurveSet<double>& curves,
                                 const std::vector<double>& times, const std::optional<Csa>& csa, const ExposureOptions& opt) {
    RL_REQUIRE(!swaps.empty(), "simulateExposure: no trades");
    // Under the textbook numeraire a fast factor makes path weights vary over many orders of magnitude;
    // measured from the last cashflow date every deflated bond stays well scaled (models/g2pp.hpp).
    const double horizon = horizonOf(swaps, curves);
    Dynamics d;
    d.twoFactor = true;
    d.loading = [&](double T) {
        const models::G2pp::Loading l = model.loading(T, horizon);
        return Load{l.h1, l.h2};
    };
    d.covariance = [&](double t) {
        const models::G2pp::Covariance c = model.covariance(t);
        return Cov{c.c11, c.c12, c.c22};
    };
    return simulate(d, swaps, curves, times, csa, opt);
}

} // namespace rl::xva