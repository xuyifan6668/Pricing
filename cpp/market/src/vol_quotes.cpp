// Swaption vol quotes: validation, JSON, hashing.
#include "rl/market/vol_quotes.hpp"

#include "rl/core/error.hpp"
#include "rl/market/hash.hpp"
#include "rl/market/swaption_cube.hpp" // cubeAxisYears, cubeNodeLabel

#include <nlohmann/json.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <set>
#include <sstream>

namespace rl::market {

std::string SwaptionVolNodeQuote::id() const { return cubeNodeLabel(expiry, tenor); }

namespace {
constexpr double kMaxVol = 0.10; // 1000 bp normal: beyond this a quote is a unit error, not a market

void addAxisLabel(std::vector<core::Tenor>& axis, core::Tenor t) {
    if (std::find(axis.begin(), axis.end(), t) == axis.end()) axis.push_back(t);
}

void sortAxis(std::vector<core::Tenor>& axis) {
    std::sort(axis.begin(), axis.end(), [](core::Tenor a, core::Tenor b) { return cubeAxisYears(a) < cubeAxisYears(b); });
    for (std::size_t i = 1; i < axis.size(); ++i)
        RL_REQUIRE(cubeAxisYears(axis[i]) > cubeAxisYears(axis[i - 1]), "vol quotes: axis labels " + axis[i - 1].str() + " and " + axis[i].str() + " are the same point");
}
} // namespace

SwaptionVolQuotes::SwaptionVolQuotes(core::Date asOf, core::Currency currency, std::string index, std::vector<SwaptionVolNodeQuote> nodes)
    : asOf_(asOf), currency_(currency), index_(std::move(index)), nodes_(std::move(nodes)) {
    RL_REQUIRE(!index_.empty(), "vol quotes: empty index name");
    RL_REQUIRE(!nodes_.empty(), "vol quotes: no nodes");
    std::set<std::string> seen;
    for (const auto& n : nodes_) {
        withContext("vol node " + n.id(), [&] {
            RL_REQUIRE(cubeAxisYears(n.expiry) > 0.0 && cubeAxisYears(n.tenor) > 0.0, "expiry and tenor must be positive");
            RL_REQUIRE(seen.insert(n.id()).second, "duplicate node");
            RL_REQUIRE(std::isfinite(n.forward) && std::fabs(n.forward) < 1.0, "forward must be a finite decimal rate");
            RL_REQUIRE(std::isfinite(n.atmNormalVol) && n.atmNormalVol > 0.0 && n.atmNormalVol < kMaxVol, "ATM normal vol must lie in (0, 1000 bp)");
            for (std::size_t k = 0; k < n.smile.size(); ++k) {
                const auto& p = n.smile[k];
                RL_REQUIRE(std::isfinite(p.offset) && p.offset != 0.0, "smile offsets must be finite and non-zero (the ATM point is atm_bp)");
                RL_REQUIRE(std::isfinite(p.normalVol) && p.normalVol > 0.0 && p.normalVol < kMaxVol, "smile vols must lie in (0, 1000 bp)");
                if (k > 0) RL_REQUIRE(p.offset > n.smile[k - 1].offset, "smile offsets must be strictly increasing");
            }
        });
        addAxisLabel(expiries_, n.expiry);
        addAxisLabel(tenors_, n.tenor);
    }
    sortAxis(expiries_);
    sortAxis(tenors_);
    for (const auto& e : expiries_)
        for (const auto& t : tenors_) RL_REQUIRE(seen.count(cubeNodeLabel(e, t)) == 1, "vol quotes: grid has no node " + cubeNodeLabel(e, t));
}

const SwaptionVolNodeQuote& SwaptionVolQuotes::find(core::Tenor expiry, core::Tenor tenor) const {
    for (const auto& n : nodes_)
        if (n.expiry == expiry && n.tenor == tenor) return n;
    throw InputError("vol quotes: no node " + cubeNodeLabel(expiry, tenor));
}

void SwaptionVolQuotes::mixInto(Hasher& h) const noexcept {
    h.mix(asOf_).mix(core::toString(currency_)).mix(std::string_view(index_));
    h.mix(static_cast<std::uint64_t>(nodes_.size()));
    for (const auto& n : nodes_) {
        h.mix(std::string_view(n.expiry.str())).mix(std::string_view(n.tenor.str()));
        h.mix(n.forward).mix(n.atmNormalVol);
        h.mix(static_cast<std::uint64_t>(n.smile.size()));
        for (const auto& p : n.smile) h.mix(p.offset).mix(p.normalVol);
    }
}

SwaptionVolQuotes SwaptionVolQuotes::fromJson(std::string_view text) {
    nlohmann::json j;
    try {
        j = nlohmann::json::parse(text);
    } catch (const nlohmann::json::exception& e) {
        throw InputError(std::string("vol quotes JSON: ") + e.what());
    }
    RL_REQUIRE(j.is_object() && j.contains("asOf") && j.contains("index") && j.contains("nodes"), "vol quotes JSON: need asOf, index and nodes");
    RL_REQUIRE(j.value("unit", std::string("normal_bp")) == "normal_bp", "vol quotes JSON: only unit normal_bp is supported in v1");
    const core::Date asOf = core::Date::fromIso(j.at("asOf").get<std::string>());
    const core::Currency ccy = j.contains("currency") ? core::currencyFromString(j.at("currency").get<std::string>()) : core::Currency::USD;
    std::vector<SwaptionVolNodeQuote> nodes;
    for (const auto& e : j.at("nodes")) {
        RL_REQUIRE(e.contains("expiry") && e.contains("tenor") && e.contains("forward_pct") && e.contains("atm_bp"),
                   "vol quotes JSON: every node needs expiry, tenor, forward_pct and atm_bp");
        SwaptionVolNodeQuote n;
        n.expiry = core::Tenor::parse(e.at("expiry").get<std::string>());
        n.tenor = core::Tenor::parse(e.at("tenor").get<std::string>());
        withContext("vol node " + n.id(), [&] {
            n.forward = e.at("forward_pct").get<double>() / 100.0;
            n.atmNormalVol = e.at("atm_bp").get<double>() / 1e4;
            if (e.contains("smile")) {
                for (const auto& p : e.at("smile")) {
                    RL_REQUIRE(p.contains("offset_bp") && p.contains("vol_bp"), "smile points need offset_bp and vol_bp");
                    n.smile.push_back(SmilePointQuote{p.at("offset_bp").get<double>() / 1e4, p.at("vol_bp").get<double>() / 1e4});
                }
            }
        });
        nodes.push_back(std::move(n));
    }
    return SwaptionVolQuotes(asOf, ccy, j.at("index").get<std::string>(), std::move(nodes));
}

SwaptionVolQuotes SwaptionVolQuotes::fromFile(const std::filesystem::path& file) {
    std::ifstream in(file);
    if (!in) throw InputError("cannot open vol quotes file " + file.string());
    std::stringstream ss;
    ss << in.rdbuf();
    return withContext("vol quotes " + file.filename().string(), [&] { return fromJson(ss.str()); });
}

std::string SwaptionVolQuotes::toJson() const {
    nlohmann::ordered_json j;
    j["asOf"] = asOf_.toIso();
    j["currency"] = std::string(core::toString(currency_));
    j["index"] = index_;
    j["unit"] = "normal_bp";
    auto arr = nlohmann::ordered_json::array();
    for (const auto& n : nodes_) {
        nlohmann::ordered_json e;
        e["expiry"] = n.expiry.str();
        e["tenor"] = n.tenor.str();
        e["forward_pct"] = n.forward * 100.0;
        e["atm_bp"] = n.atmNormalVol * 1e4;
        auto smile = nlohmann::ordered_json::array();
        for (const auto& p : n.smile) smile.push_back({{"offset_bp", p.offset * 1e4}, {"vol_bp", p.normalVol * 1e4}});
        e["smile"] = std::move(smile);
        arr.push_back(std::move(e));
    }
    j["nodes"] = std::move(arr);
    return j.dump(2);
}

} // namespace rl::market
