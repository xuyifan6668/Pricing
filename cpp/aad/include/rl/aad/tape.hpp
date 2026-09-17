#pragma once
// The AAD tape: a flat record of every operation performed on active Reals, in evaluation order.
//
// Layout (structure of arrays, one thread-local tape per thread — D10 says risk scenarios never
// share a tape):
//   node i has arguments in [argStart_[i], argStart_[i+1]); for each argument j:
//     argIndex_[j]   parent node index
//     argPartial_[j] ∂(node i)/∂(parent j), the local derivative recorded at evaluation time
//   adjoints_[i]     ∂(output)/∂(node i), filled by propagate()
// Leaves (inputs) have no arguments. Constants are never on the tape (Real::isActive() == false).
//
// Memory: 12 bytes per argument + 12 per node. A 500×400 PDE with ~10 ops per grid point per step
// is ~2e6 nodes ≈ 50 MB — fine per pricing; reset() between pricings, never between scenarios on the
// same thread without rewinding (see Position).

#include <cstddef>
#include <cstdint>
#include <vector>

namespace rl::aad {

class Tape {
public:
    static constexpr std::uint32_t kPassive = 0xFFFFFFFFu; ///< "no node": a constant.

    /// The calling thread's tape. Each thread has exactly one; it is created on first use.
    static Tape& current() noexcept;

    Tape() = default;
    Tape(const Tape&) = delete;
    Tape& operator=(const Tape&) = delete;

    // --- recording (called by Real's operators; hot) ------------------------------------------
    std::uint32_t newLeaf();
    std::uint32_t record1(std::uint32_t p, double d);
    std::uint32_t record2(std::uint32_t p1, double d1, std::uint32_t p2, double d2);

    // --- inspection ---------------------------------------------------------------------------
    std::size_t size() const noexcept { return argStart_.size() - 1; }         ///< number of nodes
    std::size_t argumentCount() const noexcept { return argIndex_.size(); }
    std::size_t arity(std::uint32_t node) const noexcept { return argStart_[node + 1] - argStart_[node]; }
    std::uint32_t parent(std::uint32_t node, std::size_t k) const noexcept { return argIndex_[argStart_[node] + k]; }
    double partial(std::uint32_t node, std::size_t k) const noexcept { return argPartial_[argStart_[node] + k]; }
    bool empty() const noexcept { return size() == 0; }

    /// Drop everything. Every Real recorded so far becomes dangling; only call between pricings.
    void reset() noexcept;

    /// A bookmark: rewindTo() truncates the tape back to it (used to discard a scenario's nodes
    /// while keeping the snapshot's leaves).
    struct Position {
        std::size_t nodes = 0;
        std::size_t args = 0;
    };
    Position mark() const noexcept { return {size(), argumentCount()}; }
    void rewindTo(Position p) noexcept;

    // --- adjoints -------------------------------------------------------------------
    /// Reverse sweep: for i from last to first, for each argument j of i: adjoint[parent] += partial·adjoint[i].
    void propagate();
    /// Sweep only over nodes in [from, to] (checkpoint-free pathwise reset: MC uses this per path).
    void propagate(std::uint32_t from, std::uint32_t to);
    double adjoint(std::uint32_t node) const noexcept { return adjoints_[node]; }
    void seedAdjoint(std::uint32_t node, double value = 1.0);
    void clearAdjoints() noexcept;
    /// Zero the adjoints of nodes ≥ first only (a nested recording must not disturb an outer one).
    void clearAdjointsFrom(std::uint32_t first) noexcept;

private:
    std::vector<std::uint32_t> argStart_{0};
    std::vector<std::uint32_t> argIndex_;
    std::vector<double> argPartial_;
    std::vector<double> adjoints_;
};

} // namespace rl::aad
