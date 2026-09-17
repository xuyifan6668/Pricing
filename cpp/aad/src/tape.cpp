#include "rl/aad/tape.hpp"

#include "rl/core/error.hpp"

namespace rl::aad {

Tape& Tape::current() noexcept {
    thread_local Tape tape;
    return tape;
}

std::uint32_t Tape::newLeaf() {
    const auto id = static_cast<std::uint32_t>(size());
    argStart_.push_back(static_cast<std::uint32_t>(argIndex_.size()));
    adjoints_.push_back(0.0);
    return id;
}

std::uint32_t Tape::record1(std::uint32_t p, double d) {
    const auto id = static_cast<std::uint32_t>(size());
    if (p != kPassive) {
        argIndex_.push_back(p);
        argPartial_.push_back(d);
    }
    argStart_.push_back(static_cast<std::uint32_t>(argIndex_.size()));
    adjoints_.push_back(0.0);
    return id;
}

std::uint32_t Tape::record2(std::uint32_t p1, double d1, std::uint32_t p2, double d2) {
    const auto id = static_cast<std::uint32_t>(size());
    if (p1 != kPassive) {
        argIndex_.push_back(p1);
        argPartial_.push_back(d1);
    }
    if (p2 != kPassive) {
        argIndex_.push_back(p2);
        argPartial_.push_back(d2);
    }
    argStart_.push_back(static_cast<std::uint32_t>(argIndex_.size()));
    adjoints_.push_back(0.0);
    return id;
}

void Tape::reset() noexcept {
    argStart_.assign(1, 0);
    argIndex_.clear();
    argPartial_.clear();
    adjoints_.clear();
}

void Tape::rewindTo(Position p) noexcept {
    if (p.nodes >= size()) return;
    argStart_.resize(p.nodes + 1);
    argIndex_.resize(p.args);
    argPartial_.resize(p.args);
    adjoints_.resize(p.nodes);
}

void Tape::seedAdjoint(std::uint32_t node, double value) {
    RL_REQUIRE(node < size(), "seedAdjoint: node index out of range");
    adjoints_[node] = value;
}

void Tape::clearAdjoints() noexcept { adjoints_.assign(adjoints_.size(), 0.0); }

void Tape::clearAdjointsFrom(std::uint32_t first) noexcept {
    for (std::size_t i = first; i < adjoints_.size(); ++i) adjoints_[i] = 0.0;
}

// Reverse sweep (Griewank & Walther 2008, ch. 3). Walking the tape backwards, each node hands its
// adjoint to its parents weighted by the local partial recorded at evaluation time:
//     adjoint[parent] += partial · adjoint[node].
// One sweep gives the derivative of ONE output with respect to EVERY input, at a cost of one
// multiply-add per recorded argument: that is the whole reason AAD beats bumping when there are
// many inputs. Nodes whose adjoint is zero do not influence the seeded output and are skipped.
void Tape::propagate() {
    if (empty()) return;
    propagate(static_cast<std::uint32_t>(size() - 1), 0);
}

/// The same sweep from node `from` down to node `to` (from ≥ to).
void Tape::propagate(std::uint32_t from, std::uint32_t to) {
    RL_REQUIRE(from < size() && to <= from, "propagate: node range out of bounds");
    for (std::uint32_t i = from + 1; i-- > to;) {
        const double a = adjoints_[i];
        if (a == 0.0) continue;
        for (std::uint32_t j = argStart_[i]; j < argStart_[i + 1]; ++j) adjoints_[argIndex_[j]] += argPartial_[j] * a;
    }
}

} // namespace rl::aad
