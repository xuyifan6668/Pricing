#pragma once
// Historical index fixings: (index name, fixing date) → rate, as published. DESIGN: passive doubles: a
// realised fixing is not differentiated, so this store is not templated on the scalar and the OIS
// kernel in linear/ mixes these doubles with T discount factors explicitly.
//
// "Fixing date" is the date the rate is *for* (the business day the overnight loan starts), not the
// publication date; SOFR for Monday is published Tuesday morning. On an as-of date d, fixings for all
// business days strictly before d are known; the fixing for d itself is not.
//
// Value type, mutated only while a snapshot is being assembled; the snapshot holds it by value or behind
// shared_ptr<const> and never writes to it. Thread-safe to read.
//
// JSON: {"USD-SOFR": {"2026-08-31": 4.31, ...}, "USD-FEDFUNDS": {...}}   (percent).

#include "rl/core/date.hpp"

#include <filesystem>
#include <map>
#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace rl::market {

class Hasher;

class FixingStore {
public:
    using Series = std::map<core::Date, double>;

    FixingStore() = default;

    /// Insert or overwrite (the later call wins; reference data corrections are normal).
    void add(std::string_view index, core::Date d, double rate);
    void add(std::string_view index, const std::vector<std::pair<core::Date, double>>& series);

    std::optional<double> fixing(std::string_view index, core::Date d) const noexcept;
    /// Throws InputError "no USD-SOFR fixing for 2026-09-03" if absent.
    double require(std::string_view index, core::Date d) const;
    bool has(std::string_view index) const noexcept;
    std::optional<core::Date> firstDate(std::string_view index) const noexcept;
    std::optional<core::Date> lastDate(std::string_view index) const noexcept;
    /// Fixings with from ≤ date ≤ to, ascending.
    std::vector<std::pair<core::Date, double>> range(std::string_view index, core::Date from, core::Date to) const;
    const Series& series(std::string_view index) const; // throws if unknown index
    std::vector<std::string> indices() const;

    std::size_t size() const noexcept; // total number of fixings
    bool empty() const noexcept { return data_.empty(); }

    /// Deterministic contribution to a snapshot hash: index names, dates and values in sorted order.
    void mixInto(Hasher& h) const noexcept;

    static FixingStore fromJson(std::string_view text);
    static FixingStore fromFile(const std::filesystem::path& file);
    std::string toJson() const;

    friend bool operator==(const FixingStore& a, const FixingStore& b) noexcept { return a.data_ == b.data_; }

private:
    std::map<std::string, Series, std::less<>> data_;
};

} // namespace rl::market
