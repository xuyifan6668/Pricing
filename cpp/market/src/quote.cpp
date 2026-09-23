#include "rl/market/quote.hpp"

#include "rl/core/error.hpp"
#include "rl/market/hash.hpp"

#include <nlohmann/json.hpp>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <sstream>
#include <unordered_set>

namespace rl::market {

namespace {
std::string lower(std::string_view s) {
    std::string u;
    for (char c : s) u.push_back(static_cast<char>(std::tolower(static_cast<unsigned char>(c))));
    return u;
}
} // namespace

std::string_view toString(QuoteType t) {
    switch (t) {
    case QuoteType::Deposit: return "deposit";
    case QuoteType::Fra: return "fra";
    case QuoteType::Future: return "future";
    case QuoteType::OisSwap: return "ois_swap";
    case QuoteType::BasisSwap: return "basis_swap";
    case QuoteType::Turn: return "turn";
    }
    return "?";
}

std::string_view toString(QuoteUnit u) {
    switch (u) {
    case QuoteUnit::Decimal: return "decimal";
    case QuoteUnit::Percent: return "pct";
    case QuoteUnit::BasisPoints: return "bp";
    case QuoteUnit::Price: return "price";
    }
    return "?";
}

QuoteType quoteTypeFromString(std::string_view s) {
    const std::string u = lower(s);
    if (u == "deposit" || u == "depo" || u == "cash") return QuoteType::Deposit;
    if (u == "fra") return QuoteType::Fra;
    if (u == "future" || u == "futures" || u == "fut") return QuoteType::Future;
    if (u == "ois_swap" || u == "ois" || u == "swap") return QuoteType::OisSwap;
    if (u == "basis_swap" || u == "basis") return QuoteType::BasisSwap;
    if (u == "turn" || u == "turn_of_year") return QuoteType::Turn;
    throw InputError("unknown quote type: '" + std::string(s) + "'");
}

QuoteUnit quoteUnitFromString(std::string_view s) {
    const std::string u = lower(s);
    if (u == "decimal" || u == "rate" || u == "dec") return QuoteUnit::Decimal;
    if (u == "pct" || u == "percent" || u == "%") return QuoteUnit::Percent;
    if (u == "bp" || u == "bps" || u == "basis_points") return QuoteUnit::BasisPoints;
    if (u == "price" || u == "px") return QuoteUnit::Price;
    throw InputError("unknown quote unit: '" + std::string(s) + "'");
}

double Quote::unitScale() const noexcept {
    switch (unit) {
    case QuoteUnit::Decimal: return 1.0;
    case QuoteUnit::Percent: return 100.0;
    case QuoteUnit::BasisPoints: return 1e4;
    case QuoteUnit::Price: return 1.0;
    }
    return 1.0;
}

double Quote::asRate() const {
    if (unit == QuoteUnit::Price) throw InputError("quote " + id + " is a price, not a rate");
    return value / unitScale();
}

double Quote::asPrice() const {
    if (unit != QuoteUnit::Price) throw InputError("quote " + id + " is a rate, not a price");
    return value;
}

double Quote::fromRate(double rate) const { return unit == QuoteUnit::Price ? rate : rate * unitScale(); }

void Quote::validate() const {
    RL_REQUIRE(!id.empty(), "quote with empty id");
    RL_REQUIRE(std::isfinite(value), "quote " + id + ": non-finite value");
    RL_REQUIRE(!index.empty(), "quote " + id + ": missing index");
    switch (type) {
    case QuoteType::Deposit:
    case QuoteType::OisSwap:
        RL_REQUIRE(tenor.has_value() || (start && end), "quote " + id + ": needs a tenor or start/end dates");
        RL_REQUIRE(unit != QuoteUnit::Price, "quote " + id + ": rate quote cannot be a price");
        break;
    case QuoteType::BasisSwap:
        RL_REQUIRE(tenor.has_value() || (start && end), "quote " + id + ": needs a tenor or start/end dates");
        RL_REQUIRE(!referenceIndex.empty(), "quote " + id + ": basis swap needs reference_index");
        RL_REQUIRE(referenceIndex != index, "quote " + id + ": basis swap legs must be on different indices");
        RL_REQUIRE(unit != QuoteUnit::Price, "quote " + id + ": spread quote cannot be a price");
        break;
    case QuoteType::Fra:
        RL_REQUIRE(start && end, "quote " + id + ": FRA needs start and end");
        RL_REQUIRE(*end > *start, "quote " + id + ": end must follow start");
        RL_REQUIRE(unit != QuoteUnit::Price, "quote " + id + ": rate quote cannot be a price");
        break;
    case QuoteType::Future:
        RL_REQUIRE(start && end, "quote " + id + ": future needs start and end (IMM reference dates)");
        RL_REQUIRE(*end > *start, "quote " + id + ": end must follow start");
        RL_REQUIRE(unit == QuoteUnit::Price, "quote " + id + ": futures are quoted as prices");
        RL_REQUIRE(value > 50.0 && value < 150.0, "quote " + id + ": futures price out of range");
        RL_REQUIRE(std::isfinite(convexityAdjustmentBp), "quote " + id + ": non-finite convexity adjustment");
        break;
    case QuoteType::Turn:
        RL_REQUIRE(start && end, "quote " + id + ": turn needs start and end");
        RL_REQUIRE(*end > *start, "quote " + id + ": end must follow start");
        RL_REQUIRE(unit != QuoteUnit::Price, "quote " + id + ": turn is a rate jump, not a price");
        break;
    }
    if (unit != QuoteUnit::Price) {
        const double r = asRate();
        RL_REQUIRE(std::fabs(r) < 1.0, "quote " + id + ": rate |" + std::to_string(r) + "| ≥ 100% — wrong unit?");
    }
}

QuoteSet::QuoteSet(core::Date asOf, core::Currency currency, std::vector<Quote> quotes)
    : asOf_(asOf), currency_(currency), quotes_(std::move(quotes)) {
    std::unordered_set<std::string> seen;
    for (const auto& q : quotes_) {
        q.validate();
        RL_REQUIRE(seen.insert(q.id).second, "duplicate quote id " + q.id);
    }
}

bool QuoteSet::contains(std::string_view id) const noexcept {
    return std::any_of(quotes_.begin(), quotes_.end(), [&](const Quote& q) { return q.id == id; });
}

std::size_t QuoteSet::indexOf(std::string_view id) const {
    for (std::size_t i = 0; i < quotes_.size(); ++i)
        if (quotes_[i].id == id) return i;
    throw InputError("no quote with id " + std::string(id));
}

const Quote& QuoteSet::find(std::string_view id) const { return quotes_[indexOf(id)]; }

std::vector<Quote> QuoteSet::ofType(QuoteType t) const {
    std::vector<Quote> out;
    for (const auto& q : quotes_)
        if (q.type == t) out.push_back(q);
    return out;
}

std::vector<std::string> QuoteSet::ids() const {
    std::vector<std::string> out;
    out.reserve(quotes_.size());
    for (const auto& q : quotes_) out.push_back(q.id);
    return out;
}

QuoteSet QuoteSet::withValues(const std::vector<double>& values) const {
    RL_REQUIRE(values.size() == quotes_.size(), "withValues: size mismatch");
    std::vector<Quote> qs = quotes_;
    for (std::size_t i = 0; i < qs.size(); ++i) qs[i].value = values[i];
    return QuoteSet(asOf_, currency_, std::move(qs));
}

QuoteSet QuoteSet::withValue(std::string_view id, double value) const {
    std::vector<Quote> qs = quotes_;
    qs[indexOf(id)].value = value;
    return QuoteSet(asOf_, currency_, std::move(qs));
}

QuoteSet QuoteSet::bumped(std::string_view id, double bump) const {
    std::vector<Quote> qs = quotes_;
    qs[indexOf(id)].value += bump;
    return QuoteSet(asOf_, currency_, std::move(qs));
}

void QuoteSet::mixInto(Hasher& h) const noexcept {
    h.mix(asOf_);
    h.mix(core::toString(currency_));
    h.mix(static_cast<std::uint64_t>(quotes_.size()));
    for (const auto& q : quotes_) {
        h.mix(std::string_view(q.id));
        h.mix(std::string_view(q.curve));
        h.mix(static_cast<int>(q.type));
        h.mix(static_cast<int>(q.unit));
        h.mix(q.value);
        h.mix(q.convexityAdjustmentBp);
    }
}

// --- JSON --------------------------------------------------------------------------------------

QuoteSet QuoteSet::fromJson(std::string_view text) {
    nlohmann::json j;
    try {
        j = nlohmann::json::parse(text);
    } catch (const nlohmann::json::exception& e) {
        throw InputError(std::string("quotes JSON: ") + e.what());
    }
    RL_REQUIRE(j.is_object() && j.contains("asOf") && j.contains("quotes"), "quotes JSON: need asOf and quotes");
    const core::Date asOf = core::Date::fromIso(j.at("asOf").get<std::string>());
    const core::Currency ccy = j.contains("currency") ? core::currencyFromString(j.at("currency").get<std::string>()) : core::Currency::USD;
    RL_REQUIRE(j.at("quotes").is_array(), "quotes JSON: quotes must be an array");
    std::vector<Quote> quotes;
    for (const auto& e : j.at("quotes")) {
        Quote q;
        RL_REQUIRE(e.contains("id") && e.contains("type") && e.contains("value"), "quotes JSON: every quote needs id, type, value");
        q.id = e.at("id").get<std::string>();
        withContext("quote " + q.id, [&] {
            q.type = quoteTypeFromString(e.at("type").get<std::string>());
            q.index = e.value("index", std::string{});
            q.referenceIndex = e.value("reference_index", std::string{});
            q.curve = e.value("curve", std::string{});
            if (e.contains("tenor")) q.tenor = core::Tenor::parse(e.at("tenor").get<std::string>());
            if (e.contains("start")) q.start = core::Date::fromIso(e.at("start").get<std::string>());
            if (e.contains("end")) q.end = core::Date::fromIso(e.at("end").get<std::string>());
            RL_REQUIRE(e.at("value").is_number(), "value must be a number");
            q.value = e.at("value").get<double>();
            q.unit = e.contains("unit") ? quoteUnitFromString(e.at("unit").get<std::string>()) : QuoteUnit::Percent;
            q.convexityAdjustmentBp = e.value("convexity_bp", 0.0);
        });
        quotes.push_back(std::move(q));
    }
    return QuoteSet(asOf, ccy, std::move(quotes));
}

QuoteSet QuoteSet::fromFile(const std::filesystem::path& file) {
    std::ifstream in(file);
    if (!in) throw InputError("cannot open quotes file " + file.string());
    std::stringstream ss;
    ss << in.rdbuf();
    return withContext("quotes " + file.filename().string(), [&] { return fromJson(ss.str()); });
}

std::string QuoteSet::toJson() const {
    nlohmann::ordered_json j;
    j["asOf"] = asOf_.toIso();
    j["currency"] = std::string(core::toString(currency_));
    auto arr = nlohmann::ordered_json::array();
    for (const auto& q : quotes_) {
        nlohmann::ordered_json e;
        e["id"] = q.id;
        e["type"] = std::string(toString(q.type));
        e["index"] = q.index;
        if (!q.referenceIndex.empty()) e["reference_index"] = q.referenceIndex;
        if (!q.curve.empty()) e["curve"] = q.curve;
        if (q.tenor) e["tenor"] = q.tenor->str();
        if (q.start) e["start"] = q.start->toIso();
        if (q.end) e["end"] = q.end->toIso();
        e["value"] = q.value;
        e["unit"] = std::string(toString(q.unit));
        if (q.type == QuoteType::Future) e["convexity_bp"] = q.convexityAdjustmentBp;
        arr.push_back(std::move(e));
    }
    j["quotes"] = std::move(arr);
    return j.dump(2);
}

} // namespace rl::market
