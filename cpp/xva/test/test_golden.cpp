// The golden regression pack: one number from every layer, stored, and compared on every build.
//
// The other tests check properties (limits, identities, agreement between methods). A change can keep
// every property and still move a price: a convention, a default, a grid size, a seed. This test pins
// the numbers themselves. Each entry has its own tolerance: tight for closed forms and deterministic
// numerics, looser for Monte Carlo, whose results can differ in the last bits between compilers.
//
// The pack is tied to the input data by a hash of the quote, vol and fixing files (line endings
// normalised). If the inputs change, the test says so instead of reporting forty failures.
//
// To regenerate after an intended change:   RL_UPDATE_GOLDEN=1 ./rl_xva_tests "[golden]"
// and commit data/golden/regression.json with the change that moved the numbers.
#include "rl/core/calendar.hpp"
#include "rl/market/hash.hpp"
#include "rl/models/g2pp_calibration.hpp"
#include "rl/models/g2pp_swaption.hpp"
#include "rl/models/lgm_calibration.hpp"
#include "rl/models/lgm_swaption.hpp"
#include "rl/pricers/bermudan_g2pp.hpp"
#include "rl/pricers/bermudan_mc.hpp"
#include "rl/pricers/bermudan_pde.hpp"
#include "rl/pricers/cheyette_mc.hpp"
#include "rl/pricers/cms_spread.hpp"
#include "rl/pricers/hedging.hpp"
#include "rl/vol/cube.hpp"
#include "rl/vol/sabr_mc.hpp"
#include "rl/xva/bermudan_exposure.hpp"
#include "rl/xva/cva_risk.hpp"
#include "rl/xva/exposure_engine.hpp"
#include "rl/xva/wrong_way.hpp"
#include "rl/xva/xva.hpp"
#include "test_support.hpp"

#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <map>
#include <sstream>
#include <string>

using namespace rl;
using namespace rl::xva::test;

namespace {

struct Entry {
    double value;
    double tolerance; ///< absolute, in the entry's own units
};

const double kExact = 1e-11, kNumeric = 1e-9, kMonteCarlo = 1e-7;

std::string inputsHash() {
    market::Hasher h;
    for (const char* f : {"quotes/usd_2026-09-04.json", "vol/usd_swaption_cube_2026-09-04.json", "fixings/usd_2026-09-04.json"}) {
        std::ifstream in(core::dataDir() / f, std::ios::binary);
        REQUIRE(in.good());
        std::stringstream ss;
        ss << in.rdbuf();
        std::string s = ss.str();
        std::erase(s, '\r'); // a Windows checkout may convert line endings
        h.mix(std::string_view(f)).mix(std::string_view(s));
    }
    return market::toHex(h.value());
}

std::vector<models::ExerciseRight> rightsOf(const instruments::Swap& swap) {
    std::vector<models::ExerciseRight> r;
    for (core::Date d : models::coterminalExerciseDates(swap, 1.0, 0, "USD-SIFMA")) r.push_back(models::makeExerciseRight(swap, usdMarket().curves, d));
    return r;
}

std::map<std::string, Entry> compute() {
    std::map<std::string, Entry> g;
    const models::UsdMarket& m = usdMarket();
    const auto& disc = m.curves.discount();

    // Curves
    for (double t : {1.0, 5.0, 10.0, 30.0}) g["curve.df." + std::to_string(static_cast<int>(t)) + "y"] = {disc.df(t), kExact};
    const instruments::Swap ten = spotSwap(10, 0.035);
    const auto rights = rightsOf(ten);
    g["curve.par_rate.10y_from_1y"] = {rights.front().forwardRate(), kExact};

    // Vol cube
    const double f5 = 0.0365; // a fixed forward: the query is the cube, not the curve
    for (int off : {-100, 0, 100}) g["vol.cube_5yx10y." + std::to_string(off) + "bp"] = {vol::cubeNormalVol<double>(m.cube, 5.0, 10.0, f5 + off * 1e-4, f5), kNumeric};
    const auto& node = m.cube.node(m.cube.expiryCount() - 1, m.cube.tenorCount() - 1);
    g["vol.node_10yx30y.rho"] = {node.sabr.rho, kNumeric};
    g["vol.node_10yx30y.nu"] = {node.sabr.nu, kNumeric};

    // LGM and the Bermudan
    const auto strip = models::makeStrip(rights, m.cube);
    const models::LgmCalibration cal = models::calibrateLgm(0.03, strip);
    g["lgm.zeta_last"] = {cal.model.zetaValues().back(), kExact};
    g["lgm.european_5y"] = {models::lgmSwaptionPrice(cal.model, rights[4]), kExact};
    for (double k : {0.025, 0.035, 0.045}) {
        const pricers::CalibratedBermudan b = pricers::calibrateBermudan(m, spotSwap(10, k), 1.0, 0, 0.03);
        g["bermudan.pde." + std::to_string(static_cast<int>(std::lround(k * 1e4))) + "bp"] = {pricers::priceBermudanPde(b.problem).value, kNumeric};
    }
    const pricers::CalibratedBermudan b = pricers::calibrateBermudan(m, ten, 1.0, 0, 0.03);
    const pricers::LsmcResult lsmc = pricers::priceBermudanLsmc(b.problem, {.trainingPaths = 4000, .pricingPaths = 8000});
    g["bermudan.lsmc_lower"] = {lsmc.lower, kMonteCarlo};
    g["bermudan.ab_upper"] = {pricers::priceBermudanUpperBound(b.problem, lsmc, {.outerPaths = 100, .innerPaths = 50}).upper, kMonteCarlo};
    g["models.kappa_marked_1y"] = {models::calibrateLgmMeanReversion(strip, models::makeAtmMatrix(m.asOf, m.curves, m.cube, {1, 2, 3, 4, 5, 6, 7, 8, 9}, {1})).kappa, 1e-6};

    // Two factors
    const models::G2pp shape(0.40, 0.03, 0.012, 0.009, -0.7);
    g["g2pp.swaption_5y"] = {models::g2ppSwaptionPrice(shape, rights[4]), kNumeric};
    const pricers::CalibratedBermudanG2pp b2 = pricers::calibrateBermudanG2pp(m, ten, 1.0, 0, shape);
    g["g2pp.bermudan_lsmc"] = {pricers::priceBermudanLsmcG2pp(b2.calibration.model, b2.rights, {.trainingPaths = 4000, .pricingPaths = 8000}).lower, kMonteCarlo};
    const auto pair = models::makeAtmMatrix(m.asOf, m.curves, m.cube, {5}, {2, 10});
    g["cms.spread_atm"] = {pricers::priceCmsSpread(shape, pair[1].right, pair[0].right, pair[1].right.forwardRate() - pair[0].right.forwardRate(), true, disc.df(pair[0].right.time), 201).price, kNumeric};

    // Cheyette
    const models::Cheyette chey(0.03, 7.0, {30.0}, {0.009});
    const pricers::CheyetteOptions co{.paths = 4000, .stepsPerYear = 12, .seed = 1};
    g["cheyette.european_5y"] = {pricers::priceEuropeansCheyette(chey, {rights[4]}, co)[0].value, kMonteCarlo};
    g["cheyette.bermudan"] = {pricers::priceBermudanCheyette(chey, rights, {.trainingPaths = 4000, .pricing = co}).lower, kMonteCarlo};

    // SABR model
    g["sabr_mc.10yx30y.atm_vol"] = {vol::sabrMonteCarlo(node.sabr, node.forward, 10.0, {node.forward}, {.paths = 4000, .stepsPerYear = 20, .seed = 1}).points[0].normalVol, kMonteCarlo};

    // Hedging
    g["hedging.weekly_pnl_sd"] = {pricers::simulateDeltaHedge(b.problem.model(), b.problem.rights(), {.paths = 200, .rebalancesPerYear = 12, .gridPoints = 201}).pnlStdDev, kMonteCarlo};

    // XVA
    const xva::HazardCurve hazard = counterpartyCurve();
    for (std::size_t i = 0; i < hazard.hazards().size(); ++i) g["xva.hazard." + std::to_string(i)] = {hazard.hazards()[i], kExact};
    g["xva.swap_cva_closed"] = {xva::swapCvaClosedForm(model(), ten, m.curves, hazard, 0.4), kNumeric};
    const std::vector<double> times = xva::exposureTimes({ten}, m.curves, 0.5);
    const xva::ExposureProfile p = xva::simulateExposure(model(), {ten}, m.curves, times, xva::Csa{0.0, 0.0, 10.0 / 252.0}, {.paths = 2000, .keepPaths = true});
    g["xva.engine_cva"] = {xva::creditAdjustment(p.times, p.epe, hazard, 0.4), kMonteCarlo};
    g["xva.engine_cva_collateralised"] = {xva::creditAdjustment(p.times, p.collateralisedEpe, hazard, 0.4), kMonteCarlo};
    g["xva.wwr_cva_beta_0.5"] = {xva::wrongWayCva(p, model(), hazard, 0.4, 0.5).cva, kMonteCarlo};
    const xva::BermudanExposure be = xva::bermudanExposure(b.problem, {.paths = 2000, .gridPoints = 401});
    g["xva.bermudan_grid_price"] = {be.price, kNumeric};
    g["xva.bermudan_cva"] = {xva::creditAdjustment(be.times, be.epe, hazard, 0.4), kMonteCarlo};
    double delta = 0.0;
    for (const auto& r : xva::swapCvaDeltaAdjoint(m, model(), ten, hazard, 0.4)) delta += r.value;
    g["xva.cva_delta_total"] = {delta, kNumeric};
    return g;
}

} // namespace

TEST_CASE("Golden regression pack: every layer reproduces its stored numbers", "[golden]") {
    const std::filesystem::path file = core::dataDir() / "golden" / "regression.json";
    const std::string hash = inputsHash();
    const std::map<std::string, Entry> now = compute();

    if (const char* update = std::getenv("RL_UPDATE_GOLDEN"); update && std::string(update) == "1") {
        nlohmann::ordered_json j;
        j["_comment"] = "Generated by cpp/xva/test/test_golden.cpp with RL_UPDATE_GOLDEN=1. Regenerate only for an intended change, and commit with it.";
        j["inputs_hash"] = hash;
        for (const auto& [name, e] : now) j["entries"][name] = {{"value", e.value}, {"tolerance", e.tolerance}};
        std::filesystem::create_directories(file.parent_path());
        std::ofstream(file) << j.dump(2) << "\n";
        WARN("regenerated " << file.string() << " (" << now.size() << " entries)");
        return;
    }

    std::ifstream in(file);
    REQUIRE(in.good());
    const nlohmann::json stored = nlohmann::json::parse(in);
    INFO("the input data changed (hash " << hash << ", pack built for " << stored.at("inputs_hash").get<std::string>()
                                         << "): regenerate with RL_UPDATE_GOLDEN=1 if the change is intended");
    REQUIRE(stored.at("inputs_hash").get<std::string>() == hash);
    const auto& entries = stored.at("entries");
    CHECK(entries.size() == now.size());
    for (const auto& [name, e] : now) {
        INFO(name);
        REQUIRE(entries.contains(name));
        const double ref = entries.at(name).at("value").get<double>(), tol = entries.at(name).at("tolerance").get<double>();
        INFO("now " << e.value << ", stored " << ref << ", tolerance " << tol);
        CHECK(std::fabs(e.value - ref) <= tol);
    }
}
