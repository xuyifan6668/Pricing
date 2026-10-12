#pragma once
// SwaptionCube<T>: the calibrated swaption vol cube as market data — an expiry × tenor grid of SABR smiles,
// each stored with the market ATM normal vol and the forward it was fitted at.
//
// DESIGN: the cube's *data* lives in market and its *mathematics* in vol. The snapshot must be able to
// hold the cube (one immutable object per market state, one hash), and market cannot depend on vol, so
// the container is here. The Hagan formula, the fit, the query that turns (expiry, tenor, strike, forward)
// into a vol, and the arbitrage diagnostic are vol/. This header knows what a SABR parameter is, not what
// it means.
//
// Interpolation (the only arithmetic here, because it is a property of the stored data, not of SABR):
//   * bilinear in (expiry, tenor) on the axis year fractions, flat extrapolation beyond both axes;
//   * ATM normal vol is interpolated as total variance σ²·t_e (linear in expiry is the no-calendar-arbitrage
//     choice for ATM; in the tenor direction at fixed expiry it is the same as interpolating σ²);
//   * ρ and ν are interpolated linearly; α is NOT interpolated — vol::cubeParams re-solves α from the
//     interpolated ATM vol at the caller's forward, so the cube reprices the liquid ATM matrix exactly off-grid
//     and the smile parameters are the decoration.
//   DESIGN: rejected interpolating α (α is tied to the forward it was fitted at through the ATM constraint;
//   interpolated α misprices ATM between nodes by several tenths of a bp) and interpolating vols by strike
//   (needs a strike grid per node; the three smile parameters are the compact, smooth representation).
//
// SIMPLIFICATION: axis times are tenor year fractions (months/12, days/365), not ACT/365F from dated expiries
// and swap start dates as a production cube keyed by trade dates would use. The difference is at most a few
// days on a 1M expiry and immaterial beyond.
//
// Value type, immutable after construction, thread-safe to read. Templated on the scalar; defined in
// detail/swaption_cube.ipp, instantiated for double (src/swaption_cube.cpp) and Real (src/instantiate_real.cpp).

#include "rl/core/date.hpp"
#include "rl/core/tenor.hpp"
#include "rl/market/sabr_params.hpp"
#include "rl/math/scalar.hpp"

#include <cstddef>
#include <filesystem>
#include <string>
#include <string_view>
#include <vector>

namespace rl::market {

class Hasher;

/// Year fraction of a cube axis label: months/12 for month and year tenors, days/365 for day and week tenors.
double cubeAxisYears(core::Tenor t);

/// "1Yx10Y": the node label used in reports and error context.
std::string cubeNodeLabel(core::Tenor expiry, core::Tenor tenor);

template <math::Scalar T>
struct CubeNode {
    double forward = 0.0;     ///< forward swap rate the smile was fitted at (passive; for the report and re-fits)
    T atmNormalVol = T(0.0);  ///< market ATM normal vol, decimal (0.0085 = 85 bp)
    SabrParams<T> sabr;       ///< fitted parameters; alpha is the α at `forward`
};

/// The cube read at an arbitrary (expiry, tenor) before α is re-solved (vol::cubeParams does that).
template <math::Scalar T>
struct InterpolatedNode {
    T atmNormalVol = T(0.0);
    T rho = T(0.0);
    T nu = T(0.0);
    double beta = 0.5;
    double shift = 0.0;
};

template <math::Scalar T>
class SwaptionCube {
public:
    /// `nodes` is row-major: nodes[i * tenors.size() + j] is (expiries[i], tenors[j]). Throws InputError on:
    /// fewer than two expiries or two tenors, axes not strictly increasing in years, a node count that is not
    /// the full grid, invalid SABR parameters, a node whose β/shift differ from the cube's, a non-positive ATM
    /// vol, or (β > 0) a node forward with F + shift ≤ 0.
    SwaptionCube(core::Date asOf, std::string index, std::vector<core::Tenor> expiries, std::vector<core::Tenor> tenors,
                 std::vector<CubeNode<T>> nodes, double beta, double shift);

    core::Date asOf() const noexcept { return asOf_; }
    const std::string& index() const noexcept { return index_; }
    double beta() const noexcept { return beta_; }
    double shift() const noexcept { return shift_; }

    const std::vector<core::Tenor>& expiryLabels() const noexcept { return expiryLabels_; }
    const std::vector<core::Tenor>& tenorLabels() const noexcept { return tenorLabels_; }
    const std::vector<double>& expiries() const noexcept { return expiries_; } ///< years
    const std::vector<double>& tenors() const noexcept { return tenors_; }     ///< years
    std::size_t expiryCount() const noexcept { return expiries_.size(); }
    std::size_t tenorCount() const noexcept { return tenors_.size(); }
    std::size_t size() const noexcept { return nodes_.size(); }

    const std::vector<CubeNode<T>>& nodes() const noexcept { return nodes_; }
    const CubeNode<T>& node(std::size_t i, std::size_t j) const; ///< throws InputError out of range
    /// (i, j) of a labelled node; throws InputError if the labels are not on the axes.
    std::pair<std::size_t, std::size_t> indexOf(core::Tenor expiry, core::Tenor tenor) const;

    /// Bilinear read at (expiry, tenor) in years; flat beyond the axes. Throws InputError unless both > 0.
    InterpolatedNode<T> interpolate(double expiry, double tenor) const;

    /// Same cube with one node replaced (re-validated). What a vega-by-node bump or a single re-fit uses.
    SwaptionCube withNode(std::size_t i, std::size_t j, CubeNode<T> node) const;

    /// Deterministic hash contribution: as-of, index, policy, labels, every node's values.
    void mixInto(Hasher& h) const noexcept;

private:
    void validate();

    core::Date asOf_;
    std::string index_;
    std::vector<core::Tenor> expiryLabels_;
    std::vector<core::Tenor> tenorLabels_;
    std::vector<double> expiries_;
    std::vector<double> tenors_;
    std::vector<CubeNode<T>> nodes_;
    double beta_ = 0.5;
    double shift_ = 0.0;
};

// --- JSON (double only: a serialised cube holds numbers, not tape nodes) ---------------------------------
//   {"asOf": "2026-09-04", "index": "USD-SOFR", "beta": 0.5, "shift": 0.0, "unit": "normal_bp",
//    "expiries": ["1M", ...], "tenors": ["1Y", ...],
//    "nodes": [{"expiry": "1M", "tenor": "1Y", "forward_pct": 3.27, "atm_bp": 105.7,
//               "alpha": 0.0521, "rho": -0.061, "nu": 0.712}, ...]}
// Nodes may appear in any order; every grid point must appear exactly once.
std::string toJson(const SwaptionCube<double>& cube);
SwaptionCube<double> swaptionCubeFromJson(std::string_view text);
SwaptionCube<double> swaptionCubeFromFile(const std::filesystem::path& file);

} // namespace rl::market
