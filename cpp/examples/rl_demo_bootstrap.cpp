// rl_demo_bootstrap <command> [interp] [args...]      (notebooks 04 and 05)
//
//   reprice    [linear|mc]                      CSV: id,type,curve,maturity,quote,model,residual,residual_bp,raw_residual
//   nodes      [linear|mc]                      CSV: curve,id,t,date,y,df,zero
//   forwards   [linear|mc] [years=5] [step=7]   CSV: date,t,sofr_fwd,ff_fwd   (1-business-day simple forwards, ACT/360)
//   sequential [linear|mc]                      CSV: id,type,curve,residual_bp_sequential,residual_bp_global
//   parrisk    [linear|mc] [tenorY=10] [K=0.03] CSV: id,type,quote,par_risk_ift,par_risk_bump,diff   (per bp, 1e6 notional payer)
//   jacobian   [linear|mc]                      CSV: row_id,col,col_curve,col_t,value   (dF/dz, dense)
//   hash       [linear|mc]                      one line: the snapshot hash of the solved curves
//
// Quotes and fixings are read from $RL_DATA_DIR/quotes/usd_2026-09-04.json and fixings/usd_2026-09-04.json.
#include "rl/bootstrap/bootstrap_problem.hpp"
#include "rl/bootstrap/global_bootstrap.hpp"
#include "rl/bootstrap/instrument_set.hpp"
#include "rl/bootstrap/jacobian.hpp"
#include "rl/bootstrap/pillars.hpp"
#include "rl/bootstrap/report.hpp"
#include "rl/core/calendar.hpp"
#include "rl/core/error.hpp"
#include "rl/instruments/builder.hpp"
#include "rl/linear/forward.hpp"
#include "rl/linear/leg_pv.hpp"
#include "rl/market/fixings.hpp"
#include "rl/market/hash.hpp"
#include "rl/market/quote.hpp"
#include "rl/market/snapshot_builder.hpp"

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>

using namespace rl;
using namespace rl::bootstrap;

namespace {

market::CurveInterpolation interpArg(int argc, char** argv) {
    if (argc < 3) return market::CurveInterpolation::LinearLogDf;
    return market::curveInterpolationFromString(argv[2]);
}

BootstrapProblem loadProblem(market::CurveInterpolation interp) {
    market::QuoteSet quotes = market::QuoteSet::fromFile(core::dataDir() / "quotes" / "usd_2026-09-04.json");
    market::FixingStore fixings = market::FixingStore::fromFile(core::dataDir() / "fixings" / "usd_2026-09-04.json");
    InstrumentSet instruments = buildInstrumentSet(quotes);
    PillarOptions po;
    po.interpolation = interp;
    std::vector<CurveSpec> specs = makeCurveSpecs(instruments, quotes.asOf(), po);
    return BootstrapProblem(quotes.asOf(), std::move(specs), std::move(instruments), std::move(fixings), std::move(quotes));
}

BootstrapResult solved(const BootstrapProblem& p) {
    BootstrapResult r = solve(p);
    if (!r.converged) throw NumericalError("bootstrap did not converge: " + (r.log.empty() ? std::string("(no log)") : r.log.back()));
    return r;
}

int usage() {
    std::cerr << "usage: rl_demo_bootstrap reprice|nodes|forwards|sequential|parrisk|jacobian|hash [linear|mc] [args]\n";
    return 2;
}

} // namespace

int main(int argc, char** argv) {
    try {
        if (argc < 2) return usage();
        const std::string cmd = argv[1];
        const auto interp = interpArg(argc, argv);
        BootstrapProblem p = loadProblem(interp);

        if (cmd == "reprice") {
            BootstrapResult r = solved(p);
            std::fputs(reprice(p, r.z).toCsv().c_str(), stdout);
            return 0;
        }
        if (cmd == "nodes") {
            BootstrapResult r = solved(p);
            auto curves = p.curves(r.z);
            std::printf("curve,id,t,date,y,df,zero\n");
            for (const auto& spec : p.curveSpecs()) {
                const auto& c = curves.curve(spec.name);
                for (std::size_t i = 0; i < spec.pillarTimes.size(); ++i) {
                    const double t = spec.pillarTimes[i], df = c.df(t);
                    std::printf("%s,%s,%.10f,%s,%.15g,%.15g,%.15g\n", spec.name.c_str(), spec.pillarIds[i].c_str(), t, spec.pillarDates[i].toIso().c_str(),
                                c.nodeValues()[i], df, -std::log(df) / t);
                }
            }
            return 0;
        }
        if (cmd == "forwards") {
            const int years = argc > 3 ? std::atoi(argv[3]) : 5;
            const int step = argc > 4 ? std::atoi(argv[4]) : 7;
            BootstrapResult r = solved(p);
            auto curves = p.curves(r.z);
            const core::Calendar cal = core::Calendar::get("USD-SIFMA");
            const auto& sofr = curves.projection("USD-SOFR");
            const auto& ff = curves.projection("USD-FEDFUNDS");
            std::printf("date,t,sofr_fwd,ff_fwd\n");
            for (core::Date d = cal.adjust(p.asOf(), core::BusinessDayConvention::Following); d <= p.asOf().addYears(years); d = cal.adjust(d + step, core::BusinessDayConvention::Following)) {
                const core::Date n = cal.advance(d, 1);
                std::printf("%s,%.8f,%.12f,%.12f\n", d.toIso().c_str(), sofr.time(d), linear::simpleForward(sofr, d, n, core::DayCount::Act360),
                            linear::simpleForward(ff, d, n, core::DayCount::Act360));
            }
            return 0;
        }
        if (cmd == "sequential") {
            BootstrapResult seq = solveSequential(p);
            BootstrapResult glob = solved(p);
            const math::Vector s = scaledResiduals(p, seq.z), g = scaledResiduals(p, glob.z);
            std::printf("id,type,curve,residual_bp_sequential,residual_bp_global\n");
            for (std::size_t i = 0; i < p.size(); ++i) {
                const auto& bi = p.instruments()[i];
                std::printf("%s,%s,%s,%.6e,%.6e\n", bi.id.c_str(), std::string(market::toString(bi.type)).c_str(), bi.curveName.c_str(),
                            s[static_cast<Eigen::Index>(i)], g[static_cast<Eigen::Index>(i)]);
            }
            return 0;
        }
        if (cmd == "parrisk") {
            const int tenor = argc > 3 ? std::atoi(argv[3]) : 10;
            const double k = argc > 4 ? std::atof(argv[4]) : 0.03;
            BootstrapResult r = solved(p);
            auto trade = instruments::OisSwapBuilder::usdSofr(p.asOf()).id("demo").tenor(core::Tenor::years(tenor)).fixedRate(k).notional(1e6).build();
            const auto& fixings = p.fixings();
            PvFunction pv = [&](const market::CurveSet<double>& c) { return linear::pv(trade, c, fixings); };
            const math::Vector g = zeroRiskByBump(p, r.z, pv);
            const math::Vector ift = parRisk(p, r.z, g);
            const math::Vector bump = parRiskByRebuild(p, r.z, pv, 1e-6);
            std::printf("id,type,quote,par_risk_ift,par_risk_bump,diff\n");
            for (std::size_t i = 0; i < p.size(); ++i) {
                const auto& bi = p.instruments()[i];
                const auto j = static_cast<Eigen::Index>(i);
                // Per basis point of the quote (a futures price bp is 0.01 of price).
                const double unit = bi.source.unit == market::QuoteUnit::Price ? -0.01 : 1e-4;
                std::printf("%s,%s,%.10g,%.8f,%.8f,%.3e\n", bi.id.c_str(), std::string(market::toString(bi.type)).c_str(), bi.source.value, ift[j] * unit,
                            bump[j] * unit, (ift[j] - bump[j]) * unit);
            }
            return 0;
        }
        if (cmd == "jacobian") {
            BootstrapResult r = solved(p);
            const math::Matrix J = dFdz(p, r.z);
            std::printf("row_id,col,col_curve,col_t,value\n");
            for (Eigen::Index i = 0; i < J.rows(); ++i) {
                Eigen::Index col = 0;
                for (const auto& spec : p.curveSpecs()) {
                    for (std::size_t k = 0; k < spec.parameterCount(); ++k, ++col) {
                        const double t = k < spec.pillarTimes.size() ? spec.pillarTimes[k] : -1.0;
                        std::printf("%s,%lld,%s,%.8f,%.12e\n", p.instrumentIds()[static_cast<std::size_t>(i)].c_str(), static_cast<long long>(col), spec.name.c_str(), t,
                                    J(i, col));
                    }
                }
            }
            return 0;
        }
        if (cmd == "hash") {
            BootstrapResult r = solved(p);
            auto snap = market::SnapshotBuilder<double>().asOf(p.asOf()).curves(p.curves(r.z)).fixings(p.fixings()).quotes(p.quotes()).build();
            std::printf("%s\n", market::toHex(snap->hash()).c_str());
            return 0;
        }
        return usage();
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << "\n";
        return 1;
    }
}
