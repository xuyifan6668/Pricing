#pragma once
// Definitions for SwaptionCube<T>. Included only by src/swaption_cube.cpp (double) and
// src/instantiate_real.cpp (Real).
#include "rl/core/error.hpp"
#include "rl/market/hash.hpp"
#include "rl/market/swaption_cube.hpp"
#include "rl/math/interpolation.hpp" // math::detail::segment

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <utility>

namespace rl::market {

template <math::Scalar T>
SwaptionCube<T>::SwaptionCube(core::Date asOf, std::string index, std::vector<core::Tenor> expiries, std::vector<core::Tenor> tenors,
                              std::vector<CubeNode<T>> nodes, double beta, double shift)
    : asOf_(asOf), index_(std::move(index)), expiryLabels_(std::move(expiries)), tenorLabels_(std::move(tenors)), nodes_(std::move(nodes)),
      beta_(beta), shift_(shift) {
    validate();
}

template <math::Scalar T>
void SwaptionCube<T>::validate() {
    RL_REQUIRE(!index_.empty(), "SwaptionCube: empty index name");
    RL_REQUIRE(expiryLabels_.size() >= 2 && tenorLabels_.size() >= 2, "SwaptionCube " + index_ + ": need at least two expiries and two tenors");
    RL_REQUIRE(std::isfinite(beta_) && beta_ >= 0.0 && beta_ <= 1.0, "SwaptionCube " + index_ + ": beta must lie in [0, 1]");
    RL_REQUIRE(std::isfinite(shift_), "SwaptionCube " + index_ + ": shift must be finite");
    expiries_.clear();
    tenors_.clear();
    for (const auto& e : expiryLabels_) expiries_.push_back(cubeAxisYears(e));
    for (const auto& t : tenorLabels_) tenors_.push_back(cubeAxisYears(t));
    for (std::size_t i = 1; i < expiries_.size(); ++i)
        RL_REQUIRE(expiries_[i] > expiries_[i - 1], "SwaptionCube " + index_ + ": expiries not strictly increasing at " + expiryLabels_[i].str());
    for (std::size_t j = 1; j < tenors_.size(); ++j)
        RL_REQUIRE(tenors_[j] > tenors_[j - 1], "SwaptionCube " + index_ + ": tenors not strictly increasing at " + tenorLabels_[j].str());
    const std::size_t ne = expiries_.size(), nt = tenors_.size();
    RL_REQUIRE(nodes_.size() == ne * nt, "SwaptionCube " + index_ + ": " + std::to_string(nodes_.size()) + " nodes for a " + std::to_string(ne) + " x " +
                                             std::to_string(nt) + " grid");
    for (std::size_t i = 0; i < ne; ++i) {
        for (std::size_t j = 0; j < nt; ++j) {
            const CubeNode<T>& n = nodes_[i * nt + j];
            withContext("SwaptionCube " + index_ + " node " + cubeNodeLabel(expiryLabels_[i], tenorLabels_[j]), [&] {
                validateSabrParams(n.sabr);
                RL_REQUIRE(n.sabr.beta == beta_ && n.sabr.shift == shift_, "node beta/shift differ from the cube's policy");
                RL_REQUIRE(math::isFinite(n.atmNormalVol) && math::value(n.atmNormalVol) > 0.0, "ATM normal vol must be positive");
                RL_REQUIRE(std::isfinite(n.forward), "non-finite forward");
                if (beta_ > 0.0) RL_REQUIRE(n.forward + shift_ > 0.0, "F + shift must be positive when beta > 0");
            });
        }
    }
}

template <math::Scalar T>
const CubeNode<T>& SwaptionCube<T>::node(std::size_t i, std::size_t j) const {
    RL_REQUIRE(i < expiries_.size() && j < tenors_.size(), "SwaptionCube " + index_ + ": node index out of range");
    return nodes_[i * tenors_.size() + j];
}

template <math::Scalar T>
std::pair<std::size_t, std::size_t> SwaptionCube<T>::indexOf(core::Tenor expiry, core::Tenor tenor) const {
    const auto ie = std::find(expiryLabels_.begin(), expiryLabels_.end(), expiry);
    const auto it = std::find(tenorLabels_.begin(), tenorLabels_.end(), tenor);
    RL_REQUIRE(ie != expiryLabels_.end() && it != tenorLabels_.end(), "SwaptionCube " + index_ + ": no node " + cubeNodeLabel(expiry, tenor));
    return {static_cast<std::size_t>(ie - expiryLabels_.begin()), static_cast<std::size_t>(it - tenorLabels_.begin())};
}

template <math::Scalar T>
InterpolatedNode<T> SwaptionCube<T>::interpolate(double expiry, double tenor) const {
    RL_REQUIRE(std::isfinite(expiry) && expiry > 0.0 && std::isfinite(tenor) && tenor > 0.0,
               "SwaptionCube " + index_ + ": query expiry and tenor must be positive");
    // Flat extrapolation: clamp the query onto the grid. Because the ATM vol is recovered as √(var / e_c) with
    // the *clamped* expiry, a short-expiry query below the first axis point returns the first node's vol, not
    // a vol that explodes as e → 0.
    const double ec = std::clamp(expiry, expiries_.front(), expiries_.back());
    const double tc = std::clamp(tenor, tenors_.front(), tenors_.back());
    const std::size_t i = math::detail::segment(expiries_, ec);
    const std::size_t j = math::detail::segment(tenors_, tc);
    const double we = (ec - expiries_[i]) / (expiries_[i + 1] - expiries_[i]);
    const double wt = (tc - tenors_[j]) / (tenors_[j + 1] - tenors_[j]);
    const std::size_t nt = tenors_.size();
    const CubeNode<T>& n00 = nodes_[i * nt + j];
    const CubeNode<T>& n01 = nodes_[i * nt + j + 1];
    const CubeNode<T>& n10 = nodes_[(i + 1) * nt + j];
    const CubeNode<T>& n11 = nodes_[(i + 1) * nt + j + 1];
    const double w00 = (1.0 - we) * (1.0 - wt), w01 = (1.0 - we) * wt, w10 = we * (1.0 - wt), w11 = we * wt;
    const double e0 = expiries_[i], e1 = expiries_[i + 1];

    InterpolatedNode<T> out;
    // Total ATM variance σ² t_e, bilinear; then back to a vol at the clamped expiry.
    const T var = w00 * (n00.atmNormalVol * n00.atmNormalVol * e0) + w01 * (n01.atmNormalVol * n01.atmNormalVol * e0) +
                  w10 * (n10.atmNormalVol * n10.atmNormalVol * e1) + w11 * (n11.atmNormalVol * n11.atmNormalVol * e1);
    using std::sqrt;
    out.atmNormalVol = sqrt(var / ec);
    out.rho = w00 * n00.sabr.rho + w01 * n01.sabr.rho + w10 * n10.sabr.rho + w11 * n11.sabr.rho;
    out.nu = w00 * n00.sabr.nu + w01 * n01.sabr.nu + w10 * n10.sabr.nu + w11 * n11.sabr.nu;
    out.beta = beta_;
    out.shift = shift_;
    return out;
}

template <math::Scalar T>
SwaptionCube<T> SwaptionCube<T>::withNode(std::size_t i, std::size_t j, CubeNode<T> node) const {
    RL_REQUIRE(i < expiries_.size() && j < tenors_.size(), "SwaptionCube " + index_ + ": withNode index out of range");
    std::vector<CubeNode<T>> nodes = nodes_;
    nodes[i * tenors_.size() + j] = std::move(node);
    return SwaptionCube(asOf_, index_, expiryLabels_, tenorLabels_, std::move(nodes), beta_, shift_);
}

template <math::Scalar T>
void SwaptionCube<T>::mixInto(Hasher& h) const noexcept {
    h.mix(asOf_).mix(std::string_view(index_)).mix(beta_).mix(shift_);
    h.mix(static_cast<std::uint64_t>(expiryLabels_.size()));
    for (const auto& e : expiryLabels_) h.mix(std::string_view(e.str()));
    h.mix(static_cast<std::uint64_t>(tenorLabels_.size()));
    for (const auto& t : tenorLabels_) h.mix(std::string_view(t.str()));
    for (const auto& n : nodes_) {
        h.mix(n.forward).mix(math::value(n.atmNormalVol));
        mixSabrParams(h, n.sabr);
    }
}

} // namespace rl::market

#define RL_MARKET_INSTANTIATE_CUBE(T) template class rl::market::SwaptionCube<T>;
