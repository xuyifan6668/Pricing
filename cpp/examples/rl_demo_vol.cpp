// rl_demo_vol <command> [args...]      (notebooks 06 and 07)
//
// All numeric inputs are decimals (0.035 = 3.5%, 0.0085 = 85 bp); output columns carry their units.
//
//   convert   F K T sigmaN shift                         CSV: forward,strike,expiry,shift,normal_vol_bp,lognormal_vol_pct,roundtrip_normal_vol_bp,atm_series_lognormal_pct
//   smile     alpha beta rho nu shift F T [kMin kMax n]  CSV: strike,offset_bp,normal_vol_bp,lognormal_vol_pct
//   density   alpha beta rho nu shift F T [kMin kMax n]  CSV: strike,normal_vol_bp,call,dcdk,density
//   check     alpha beta rho nu shift F T [nStdDev n]    CSV (one row): k_min,k_max,min_density,min_density_strike,first_negative_below,first_negative_above,call_spread_ok,nonpositive_vol,arbitrage_free,mass,mean
//   quotes                                               CSV: id,expiry,tenor,expiry_years,tenor_years,forward_pct,offset_bp,strike_pct,vol_bp
//   fit       all|EXPIRY TENOR [beta=0.5 shift=0]        CSV: id,expiry,tenor,forward_pct,atm_bp,alpha,beta,rho,nu,rms_bp,max_bp,iterations,status,ok,flags
//   residuals [beta shift]                               CSV: id,expiry,tenor,quote,strike_pct,market_bp,model_bp,residual_bp
//   cube      [beta shift]                               the fitted cube as JSON (market::toJson)
//   arbitrage [beta shift nStdDev=4 n=801]               CSV: id,expiry,tenor,forward_pct,atm_bp,k_min_pct,k_max_pct,min_density,first_negative_below_pct,first_negative_above_pct,call_spread_ok,nonpositive_vol,arbitrage_free
//   query     cubeFile expiry tenor F kMin kMax n        CSV: strike,normal_vol_bp   (expiry, tenor in years)
//
#include "rl/core/calendar.hpp"
#include "rl/core/error.hpp"
#include "rl/core/tenor.hpp"
#include "rl/market/swaption_cube.hpp"
#include "rl/market/vol_quotes.hpp"
#include "rl/vol/arbitrage.hpp"
#include "rl/vol/conversion.hpp"
#include "rl/vol/cube.hpp"
#include "rl/vol/sabr.hpp"
#include "rl/vol/sabr_fit.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <optional>
#include <string>

using namespace rl;

namespace {

double arg(int argc, char** argv, int i, const char* name) {
    if (argc <= i) throw InputError(std::string("missing argument: ") + name);
    return std::atof(argv[i]);
}
double argOr(int argc, char** argv, int i, double fallback) { return argc > i ? std::atof(argv[i]) : fallback; }

market::SabrParams<double> paramsArg(int argc, char** argv, int i) {
    market::SabrParams<double> p;
    p.alpha = arg(argc, argv, i, "alpha");
    p.beta = arg(argc, argv, i + 1, "beta");
    p.rho = arg(argc, argv, i + 2, "rho");
    p.nu = arg(argc, argv, i + 3, "nu");
    p.shift = arg(argc, argv, i + 4, "shift");
    market::validateSabrParams(p);
    return p;
}

vol::SabrFitOptions policyArg(int argc, char** argv, int i) {
    vol::SabrFitOptions o;
    o.beta = argOr(argc, argv, i, 0.5);
    o.shift = argOr(argc, argv, i + 1, 0.0);
    return o;
}

market::SwaptionVolQuotes loadQuotes() { return market::SwaptionVolQuotes::fromFile(core::dataDir() / "vol" / "usd_swaption_cube_2026-09-04.json"); }

std::string joinFlags(const std::vector<std::string>& flags) {
    std::string out;
    for (const auto& f : flags) {
        if (!out.empty()) out += ";";
        for (char c : f) out.push_back(c == ',' ? ' ' : c);
    }
    return out;
}

std::string pct(const std::optional<double>& v) {
    if (!v) return "";
    char buf[32];
    std::snprintf(buf, sizeof buf, "%.6f", *v * 100.0);
    return buf;
}

/// Strike range for the smile/density/check commands: F ± nStdDev·σ_ATM√T, floored for β > 0.
std::vector<double> strikesArg(int argc, char** argv, int i, const market::SabrParams<double>& p, double f, double t, double nStdDev) {
    const int n = argc > i + 2 ? std::atoi(argv[i + 2]) : 201;
    if (argc > i + 1) return vol::strikeGrid(arg(argc, argv, i, "kMin"), arg(argc, argv, i + 1, "kMax"), static_cast<std::size_t>(n));
    const auto [lo, hi] = vol::diagnosticStrikeRange(p, f, t, vol::sabrAtmNormalVol(p, f, t), nStdDev);
    return vol::strikeGrid(lo, hi, static_cast<std::size_t>(n));
}

std::vector<const market::SwaptionVolNodeQuote*> gridNodes(const market::SwaptionVolQuotes& q) {
    std::vector<const market::SwaptionVolNodeQuote*> out;
    for (const auto& e : q.expiries())
        for (const auto& t : q.tenors()) out.push_back(&q.find(e, t));
    return out;
}

int usage() {
    std::cerr << "usage: rl_demo_vol convert|smile|density|check|quotes|fit|residuals|cube|arbitrage|query [args]  (see the header of rl_demo_vol.cpp)\n";
    return 2;
}

} // namespace

int main(int argc, char** argv) {
    try {
        if (argc < 2) return usage();
        const std::string cmd = argv[1];

        if (cmd == "convert") {
            const double f = arg(argc, argv, 2, "F"), k = arg(argc, argv, 3, "K"), t = arg(argc, argv, 4, "T");
            const double sn = arg(argc, argv, 5, "sigmaN"), s = arg(argc, argv, 6, "shift");
            const double sl = vol::normalToShiftedLognormalVol(f, k, t, sn, s);
            const double back = vol::shiftedLognormalToNormalVol(f, k, t, sl, s);
            std::printf("forward,strike,expiry,shift,normal_vol_bp,lognormal_vol_pct,roundtrip_normal_vol_bp,atm_series_lognormal_pct\n");
            std::printf("%.10g,%.10g,%.10g,%.10g,%.12g,%.12g,%.12g,%.12g\n", f, k, t, s, sn * 1e4, sl * 100.0, back * 1e4,
                        vol::atmLognormalFromNormalApprox(f, t, sn, s) * 100.0);
            return 0;
        }
        if (cmd == "smile" || cmd == "density" || cmd == "check") {
            const auto p = paramsArg(argc, argv, 2);
            const double f = arg(argc, argv, 7, "F"), t = arg(argc, argv, 8, "T");
            if (cmd == "check") {
                const double nsd = argOr(argc, argv, 9, 4.0);
                const int n = argc > 10 ? std::atoi(argv[10]) : 801;
                const auto [lo, hi] = vol::diagnosticStrikeRange(p, f, t, vol::sabrAtmNormalVol(p, f, t), nsd);
                const auto c = vol::checkSabrSmile(p, f, t, lo, hi, static_cast<std::size_t>(n));
                std::printf("k_min,k_max,min_density,min_density_strike,first_negative_below,first_negative_above,call_spread_ok,nonpositive_vol,arbitrage_free,mass,mean\n");
                auto opt = [](const std::optional<double>& v) { return v ? std::to_string(*v) : std::string(); };
                std::printf("%.10g,%.10g,%.10g,%.10g,%s,%s,%d,%d,%d,%.10g,%.10g\n", c.kMin, c.kMax, c.minDensity, c.minDensityStrike,
                            opt(c.firstNegativeDensityBelowForward).c_str(), opt(c.firstNegativeDensityAboveForward).c_str(), c.callSpreadViolated ? 0 : 1,
                            c.nonPositiveVol ? 1 : 0, c.arbitrageFree() ? 1 : 0, c.integratedMass, c.impliedMean);
                return 0;
            }
            const auto ks = strikesArg(argc, argv, 9, p, f, t, 3.0);
            if (cmd == "smile") {
                std::printf("strike,offset_bp,normal_vol_bp,lognormal_vol_pct\n");
                for (double k : ks) {
                    const double vn = vol::sabrNormalVol(p, f, k, t);
                    const bool hasLn = f + p.shift > 0.0 && k + p.shift > 0.0;
                    const double vl = hasLn ? vol::sabrLognormalVol(p, f, k, t) : std::nan("");
                    std::printf("%.10g,%.6f,%.10g,%.10g\n", k, (k - f) * 1e4, vn * 1e4, vl * 100.0);
                }
                return 0;
            }
            std::printf("strike,normal_vol_bp,call,dcdk,density\n");
            for (const auto& d : vol::riskNeutralDensity(p, f, t, ks))
                std::printf("%.10g,%.10g,%.12g,%.12g,%.12g\n", d.strike, d.normalVol * 1e4, d.callPrice, d.dCalldK, d.density);
            return 0;
        }
        if (cmd == "quotes") {
            const auto q = loadQuotes();
            std::printf("id,expiry,tenor,expiry_years,tenor_years,forward_pct,offset_bp,strike_pct,vol_bp\n");
            for (const auto* n : gridNodes(q)) {
                for (const auto& sq : vol::smileQuotes(*n)) {
                    std::printf("%s,%s,%s,%.10g,%.10g,%.6f,%.4f,%.6f,%.4f\n", n->id().c_str(), n->expiry.str().c_str(), n->tenor.str().c_str(),
                                market::cubeAxisYears(n->expiry), market::cubeAxisYears(n->tenor), n->forward * 100.0, (sq.strike - n->forward) * 1e4,
                                sq.strike * 100.0, sq.normalVol * 1e4);
                }
            }
            return 0;
        }
        if (cmd == "fit") {
            const auto q = loadQuotes();
            std::vector<const market::SwaptionVolNodeQuote*> nodes;
            int policyAt = 3;
            if (argc > 2 && std::string(argv[2]) == "all") {
                nodes = gridNodes(q);
            } else {
                if (argc < 4) return usage();
                nodes.push_back(&q.find(core::Tenor::parse(argv[2]), core::Tenor::parse(argv[3])));
                policyAt = 4;
            }
            const auto opt = policyArg(argc, argv, policyAt);
            std::printf("id,expiry,tenor,forward_pct,atm_bp,alpha,beta,rho,nu,rms_bp,max_bp,iterations,status,ok,flags\n");
            for (const auto* n : nodes) {
                const auto r = withContext("node " + n->id(), [&] { return vol::fitSabrSmile(vol::makeNodeProblem(*n, opt)); });
                std::printf("%s,%s,%s,%.6f,%.4f,%.12g,%.4g,%.12g,%.12g,%.6e,%.6e,%d,%s,%d,%s\n", n->id().c_str(), n->expiry.str().c_str(), n->tenor.str().c_str(),
                            n->forward * 100.0, n->atmNormalVol * 1e4, r.params.alpha, r.params.beta, r.params.rho, r.params.nu, r.rmsBp, r.maxAbsBp, r.iterations,
                            math::toString(r.status), r.ok ? 1 : 0, joinFlags(r.flags).c_str());
            }
            return 0;
        }
        if (cmd == "residuals") {
            const auto q = loadQuotes();
            const auto opt = policyArg(argc, argv, 2);
            std::printf("id,expiry,tenor,quote,strike_pct,market_bp,model_bp,residual_bp\n");
            for (const auto* n : gridNodes(q)) {
                const auto prob = vol::makeNodeProblem(*n, opt);
                const auto r = vol::fitSabrSmile(prob);
                for (std::size_t i = 0; i < prob.quotes().size(); ++i) {
                    const auto& sq = prob.quotes()[i];
                    std::printf("%s,%s,%s,%s,%.6f,%.4f,%.6f,%.6e\n", n->id().c_str(), n->expiry.str().c_str(), n->tenor.str().c_str(), prob.instrumentIds()[i].c_str(),
                                sq.strike * 100.0, sq.normalVol * 1e4, sq.normalVol * 1e4 + r.residualsBp[i], r.residualsBp[i]);
                }
            }
            return 0;
        }
        if (cmd == "cube") {
            const auto fit = vol::fitCube(loadQuotes(), policyArg(argc, argv, 2));
            std::fputs(market::toJson(fit.cube).c_str(), stdout);
            std::fputs("\n", stdout);
            return 0;
        }
        if (cmd == "arbitrage") {
            const auto q = loadQuotes();
            const auto opt = policyArg(argc, argv, 2);
            const double nsd = argOr(argc, argv, 4, 4.0);
            const int npts = argc > 5 ? std::atoi(argv[5]) : 801;
            std::printf("id,expiry,tenor,forward_pct,atm_bp,k_min_pct,k_max_pct,min_density,first_negative_below_pct,first_negative_above_pct,call_spread_ok,nonpositive_vol,arbitrage_free\n");
            for (const auto* n : gridNodes(q)) {
                const auto r = vol::fitSabrSmile(vol::makeNodeProblem(*n, opt));
                const double t = market::cubeAxisYears(n->expiry);
                const auto [lo, hi] = vol::diagnosticStrikeRange(r.params, n->forward, t, n->atmNormalVol, nsd);
                const auto c = vol::checkSabrSmile(r.params, n->forward, t, lo, hi, static_cast<std::size_t>(npts));
                std::printf("%s,%s,%s,%.6f,%.4f,%.6f,%.6f,%.6e,%s,%s,%d,%d,%d\n", n->id().c_str(), n->expiry.str().c_str(), n->tenor.str().c_str(), n->forward * 100.0,
                            n->atmNormalVol * 1e4, lo * 100.0, hi * 100.0, c.minDensity, pct(c.firstNegativeDensityBelowForward).c_str(),
                            pct(c.firstNegativeDensityAboveForward).c_str(), c.callSpreadViolated ? 0 : 1, c.nonPositiveVol ? 1 : 0, c.arbitrageFree() ? 1 : 0);
            }
            return 0;
        }
        if (cmd == "query") {
            if (argc < 9) return usage();
            const auto cube = market::swaptionCubeFromFile(argv[2]);
            const double e = arg(argc, argv, 3, "expiry"), t = arg(argc, argv, 4, "tenor"), f = arg(argc, argv, 5, "F");
            const auto ks = vol::strikeGrid(arg(argc, argv, 6, "kMin"), arg(argc, argv, 7, "kMax"), static_cast<std::size_t>(std::atoi(argv[8])));
            const auto vols = vol::cubeSmile(cube, e, t, f, ks);
            std::printf("strike,normal_vol_bp\n");
            for (std::size_t i = 0; i < ks.size(); ++i) std::printf("%.10g,%.10g\n", ks[i], vols[i] * 1e4);
            return 0;
        }
        return usage();
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << "\n";
        return 1;
    }
}
