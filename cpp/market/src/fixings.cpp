#include "rl/market/fixings.hpp"

#include "rl/core/error.hpp"
#include "rl/market/hash.hpp"

#include <nlohmann/json.hpp>

#include <cmath>
#include <fstream>
#include <sstream>

namespace rl::market {

void FixingStore::add(std::string_view index, core::Date d, double rate) {
    RL_REQUIRE(!index.empty(), "fixing: empty index name");
    RL_REQUIRE(std::isfinite(rate), "fixing: non-finite rate for " + std::string(index) + " " + d.toIso());
    auto it = data_.find(index);
    if (it == data_.end()) it = data_.emplace(std::string(index), Series{}).first;
    it->second[d] = rate;
}

void FixingStore::add(std::string_view index, const std::vector<std::pair<core::Date, double>>& series) {
    for (const auto& [d, r] : series) add(index, d, r);
}

std::optional<double> FixingStore::fixing(std::string_view index, core::Date d) const noexcept {
    auto it = data_.find(index);
    if (it == data_.end()) return std::nullopt;
    auto jt = it->second.find(d);
    if (jt == it->second.end()) return std::nullopt;
    return jt->second;
}

double FixingStore::require(std::string_view index, core::Date d) const {
    auto f = fixing(index, d);
    if (!f) throw InputError("no " + std::string(index) + " fixing for " + d.toIso());
    return *f;
}

bool FixingStore::has(std::string_view index) const noexcept { return data_.find(index) != data_.end(); }

std::optional<core::Date> FixingStore::firstDate(std::string_view index) const noexcept {
    auto it = data_.find(index);
    if (it == data_.end() || it->second.empty()) return std::nullopt;
    return it->second.begin()->first;
}

std::optional<core::Date> FixingStore::lastDate(std::string_view index) const noexcept {
    auto it = data_.find(index);
    if (it == data_.end() || it->second.empty()) return std::nullopt;
    return it->second.rbegin()->first;
}

std::vector<std::pair<core::Date, double>> FixingStore::range(std::string_view index, core::Date from, core::Date to) const {
    std::vector<std::pair<core::Date, double>> out;
    auto it = data_.find(index);
    if (it == data_.end()) return out;
    for (auto jt = it->second.lower_bound(from); jt != it->second.end() && jt->first <= to; ++jt) out.emplace_back(jt->first, jt->second);
    return out;
}

const FixingStore::Series& FixingStore::series(std::string_view index) const {
    auto it = data_.find(index);
    if (it == data_.end()) throw InputError("no fixings for index " + std::string(index));
    return it->second;
}

std::vector<std::string> FixingStore::indices() const {
    std::vector<std::string> out;
    for (const auto& [k, v] : data_) out.push_back(k);
    return out;
}

std::size_t FixingStore::size() const noexcept {
    std::size_t n = 0;
    for (const auto& [k, v] : data_) n += v.size();
    return n;
}

void FixingStore::mixInto(Hasher& h) const noexcept {
    h.mix(static_cast<std::uint64_t>(data_.size()));
    for (const auto& [name, series] : data_) {
        h.mix(std::string_view(name));
        h.mix(static_cast<std::uint64_t>(series.size()));
        for (const auto& [d, r] : series) {
            h.mix(d);
            h.mix(r);
        }
    }
}

FixingStore FixingStore::fromJson(std::string_view text) {
    nlohmann::json j;
    try {
        j = nlohmann::json::parse(text);
    } catch (const nlohmann::json::exception& e) {
        throw InputError(std::string("fixings JSON: ") + e.what());
    }
    RL_REQUIRE(j.is_object(), "fixings JSON: top level must be an object keyed by index name");
    FixingStore fs;
    for (const auto& [index, series] : j.items()) {
        if (!index.empty() && index[0] == '_') continue; // "_comment" and friends
        RL_REQUIRE(series.is_object(), "fixings JSON: series for " + index + " must be an object keyed by date");
        for (const auto& [date, v] : series.items()) {
            RL_REQUIRE(v.is_number(), "fixings JSON: non-numeric fixing for " + index + " " + date);
            fs.add(index, core::Date::fromIso(date), v.get<double>() / 100.0);
        }
    }
    return fs;
}

FixingStore FixingStore::fromFile(const std::filesystem::path& file) {
    std::ifstream in(file);
    if (!in) throw InputError("cannot open fixings file " + file.string());
    std::stringstream ss;
    ss << in.rdbuf();
    return withContext("fixings " + file.filename().string(), [&] { return fromJson(ss.str()); });
}

std::string FixingStore::toJson() const {
    nlohmann::ordered_json j = nlohmann::ordered_json::object();
    for (const auto& [name, series] : data_) {
        nlohmann::ordered_json s = nlohmann::ordered_json::object();
        for (const auto& [d, r] : series) s[d.toIso()] = r * 100.0;
        j[name] = std::move(s);
    }
    return j.dump(2);
}

} // namespace rl::market
