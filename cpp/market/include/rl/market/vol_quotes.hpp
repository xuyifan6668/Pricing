#pragma once
// SwaptionVolQuotes: the swaption vol screen as it arrives — for each (expiry, tenor) node an ATM normal vol
// and a smile of normal vols at strike offsets from the forward. No SABR, no fit; vol/sabr_fit.hpp and the
// Python cube builder turn this into a SwaptionCube (docs/DESIGN.md §3: cube *assembly* is Python, the
// per-node fit kernel is C++).
//
// SIMPLIFICATION: each node carries its forward swap rate. Production computes the forward from the curve
// at build time (linear::parRate on the node's underlying swap), so a vol screen and a curve built from
// different snapshots are caught. Storing it keeps testable without a solved bootstrap and makes a
// node fit reproducible from one file. the calibration recomputes forwards from the curve and reads
// only vols from the cube.
//
// DESIGN: no β, shift or weighting in the quote file. Those are model policy, chosen by the cube builder
// (python/ratelib/vol/cube_builder.py) and passed to the fit; a quote file that carried them would make a
// policy change look like a market move in the snapshot hash.
//
// JSON (data/vol/usd_swaption_cube_2026-09-04.json), human units:
//   {"asOf": "2026-09-04", "currency": "USD", "index": "USD-SOFR", "unit": "normal_bp",
//    "nodes": [{"expiry": "1Y", "tenor": "10Y", "forward_pct": 3.60, "atm_bp": 88.0,
//               "smile": [{"offset_bp": -100, "vol_bp": 93.1}, ..., {"offset_bp": 100, "vol_bp": 86.0}]}, ...]}
// In memory everything is decimal: forward 0.036, vols 0.0088, offsets −0.01.

#include "rl/core/currency.hpp"
#include "rl/core/date.hpp"
#include "rl/core/tenor.hpp"

#include <cstddef>
#include <filesystem>
#include <string>
#include <string_view>
#include <vector>

namespace rl::market {

class Hasher;

struct SmilePointQuote {
    double offset = 0.0;    ///< strike − forward, decimal; never 0 (the ATM point is the node's atmNormalVol)
    double normalVol = 0.0; ///< decimal
};

struct SwaptionVolNodeQuote {
    core::Tenor expiry;
    core::Tenor tenor;
    double forward = 0.0;      ///< decimal (see SIMPLIFICATION above)
    double atmNormalVol = 0.0; ///< decimal
    std::vector<SmilePointQuote> smile; ///< offsets strictly increasing

    std::string id() const;    ///< "1Yx10Y"
};

class SwaptionVolQuotes {
public:
    SwaptionVolQuotes() = default;
    /// Throws InputError on: no nodes, a duplicate (expiry, tenor), a grid with holes (the cube needs every
    /// expiry × tenor), a non-positive or implausible (≥ 1000 bp) vol, a non-finite forward, smile offsets that
    /// are not strictly increasing or contain zero. Errors carry the node id.
    SwaptionVolQuotes(core::Date asOf, core::Currency currency, std::string index, std::vector<SwaptionVolNodeQuote> nodes);

    static SwaptionVolQuotes fromJson(std::string_view text);
    static SwaptionVolQuotes fromFile(const std::filesystem::path& file);
    std::string toJson() const;

    core::Date asOf() const noexcept { return asOf_; }
    core::Currency currency() const noexcept { return currency_; }
    const std::string& index() const noexcept { return index_; }
    const std::vector<SwaptionVolNodeQuote>& nodes() const noexcept { return nodes_; }
    std::size_t size() const noexcept { return nodes_.size(); }

    /// Axis labels in increasing year fraction (cubeAxisYears).
    const std::vector<core::Tenor>& expiries() const noexcept { return expiries_; }
    const std::vector<core::Tenor>& tenors() const noexcept { return tenors_; }
    /// Throws InputError if absent.
    const SwaptionVolNodeQuote& find(core::Tenor expiry, core::Tenor tenor) const;

    /// Deterministic contribution to a snapshot hash: as-of, currency, index, every node's numbers.
    void mixInto(Hasher& h) const noexcept;

private:
    core::Date asOf_;
    core::Currency currency_ = core::Currency::USD;
    std::string index_;
    std::vector<SwaptionVolNodeQuote> nodes_;
    std::vector<core::Tenor> expiries_;
    std::vector<core::Tenor> tenors_;
};

} // namespace rl::market
