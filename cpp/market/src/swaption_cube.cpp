// SwaptionCube<double> instantiation, axis helpers and JSON.
#include "rl/market/detail/swaption_cube.ipp"

#include <nlohmann/json.hpp>

#include <fstream>
#include <optional>
#include <sstream>

namespace rl::market {

double cubeAxisYears(core::Tenor t) {
    RL_REQUIRE(t.n > 0, "cube axis tenor must be positive, got " + t.str());
    switch (t.unit) {
    case core::TenorUnit::Days: return t.n / 365.0;
    case core::TenorUnit::Weeks: return 7.0 * t.n / 365.0;
    case core::TenorUnit::Months: return t.n / 12.0;
    case core::TenorUnit::Years: return static_cast<double>(t.n);
    }
    return 0.0;
}

std::string cubeNodeLabel(core::Tenor expiry, core::Tenor tenor) { return expiry.str() + "x" + tenor.str(); }

std::string toJson(const SwaptionCube<double>& cube) {
    nlohmann::ordered_json j;
    j["asOf"] = cube.asOf().toIso();
    j["index"] = cube.index();
    j["beta"] = cube.beta();
    j["shift"] = cube.shift();
    j["unit"] = "normal_bp";
    auto ex = nlohmann::ordered_json::array();
    for (const auto& e : cube.expiryLabels()) ex.push_back(e.str());
    auto tn = nlohmann::ordered_json::array();
    for (const auto& t : cube.tenorLabels()) tn.push_back(t.str());
    j["expiries"] = std::move(ex);
    j["tenors"] = std::move(tn);
    auto nodes = nlohmann::ordered_json::array();
    for (std::size_t i = 0; i < cube.expiryCount(); ++i) {
        for (std::size_t k = 0; k < cube.tenorCount(); ++k) {
            const auto& n = cube.node(i, k);
            nlohmann::ordered_json e;
            e["expiry"] = cube.expiryLabels()[i].str();
            e["tenor"] = cube.tenorLabels()[k].str();
            e["forward_pct"] = n.forward * 100.0;
            e["atm_bp"] = n.atmNormalVol * 1e4;
            e["alpha"] = n.sabr.alpha;
            e["rho"] = n.sabr.rho;
            e["nu"] = n.sabr.nu;
            nodes.push_back(std::move(e));
        }
    }
    j["nodes"] = std::move(nodes);
    return j.dump(2);
}

SwaptionCube<double> swaptionCubeFromJson(std::string_view text) {
    nlohmann::json j;
    try {
        j = nlohmann::json::parse(text);
    } catch (const nlohmann::json::exception& e) {
        throw InputError(std::string("swaption cube JSON: ") + e.what());
    }
    RL_REQUIRE(j.is_object(), "swaption cube JSON: top level must be an object");
    for (const char* key : {"asOf", "index", "beta", "expiries", "tenors", "nodes"})
        RL_REQUIRE(j.contains(key), std::string("swaption cube JSON: missing '") + key + "'");
    RL_REQUIRE(j.value("unit", std::string("normal_bp")) == "normal_bp", "swaption cube JSON: only unit normal_bp is supported");
    const core::Date asOf = core::Date::fromIso(j.at("asOf").get<std::string>());
    const std::string index = j.at("index").get<std::string>();
    const double beta = j.at("beta").get<double>();
    const double shift = j.value("shift", 0.0);
    std::vector<core::Tenor> ex, tn;
    for (const auto& e : j.at("expiries")) ex.push_back(core::Tenor::parse(e.get<std::string>()));
    for (const auto& t : j.at("tenors")) tn.push_back(core::Tenor::parse(t.get<std::string>()));
    RL_REQUIRE(!ex.empty() && !tn.empty(), "swaption cube JSON: empty axes");

    std::vector<std::optional<CubeNode<double>>> slots(ex.size() * tn.size());
    for (const auto& e : j.at("nodes")) {
        RL_REQUIRE(e.contains("expiry") && e.contains("tenor"), "swaption cube JSON: every node needs expiry and tenor");
        const core::Tenor te = core::Tenor::parse(e.at("expiry").get<std::string>());
        const core::Tenor tt = core::Tenor::parse(e.at("tenor").get<std::string>());
        const std::string label = cubeNodeLabel(te, tt);
        withContext("node " + label, [&] {
            const auto ie = std::find(ex.begin(), ex.end(), te);
            const auto it = std::find(tn.begin(), tn.end(), tt);
            RL_REQUIRE(ie != ex.end() && it != tn.end(), "node is not on the declared axes");
            auto& slot = slots[static_cast<std::size_t>(ie - ex.begin()) * tn.size() + static_cast<std::size_t>(it - tn.begin())];
            RL_REQUIRE(!slot.has_value(), "duplicate node");
            for (const char* key : {"forward_pct", "atm_bp", "alpha", "rho", "nu"}) RL_REQUIRE(e.contains(key), std::string("missing '") + key + "'");
            CubeNode<double> n;
            n.forward = e.at("forward_pct").get<double>() / 100.0;
            n.atmNormalVol = e.at("atm_bp").get<double>() / 1e4;
            n.sabr.alpha = e.at("alpha").get<double>();
            n.sabr.rho = e.at("rho").get<double>();
            n.sabr.nu = e.at("nu").get<double>();
            n.sabr.beta = beta;
            n.sabr.shift = shift;
            slot = n;
        });
    }
    std::vector<CubeNode<double>> nodes;
    nodes.reserve(slots.size());
    for (std::size_t i = 0; i < ex.size(); ++i)
        for (std::size_t k = 0; k < tn.size(); ++k) {
            const auto& s = slots[i * tn.size() + k];
            RL_REQUIRE(s.has_value(), "swaption cube JSON: missing node " + cubeNodeLabel(ex[i], tn[k]));
            nodes.push_back(*s);
        }
    return SwaptionCube<double>(asOf, index, std::move(ex), std::move(tn), std::move(nodes), beta, shift);
}

SwaptionCube<double> swaptionCubeFromFile(const std::filesystem::path& file) {
    std::ifstream in(file);
    if (!in) throw InputError("cannot open swaption cube file " + file.string());
    std::stringstream ss;
    ss << in.rdbuf();
    return withContext("swaption cube " + file.filename().string(), [&] { return swaptionCubeFromJson(ss.str()); });
}

} // namespace rl::market

RL_MARKET_INSTANTIATE_CUBE(double)
