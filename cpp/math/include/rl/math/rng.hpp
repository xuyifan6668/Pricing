#pragma once
// Pseudo-random numbers and the seed policy.
//
// The engine is std::mt19937_64: its output sequence is specified by the C++ standard bit for bit,
// so MSVC and libstdc++ agree. std::uniform_real_distribution and std::normal_distribution are NOT
// specified (each library implements them differently, and libstdc++'s normal is a Box–Muller with
// a cached second value that reorders the stream), which is why this class does the uniform and the
// gaussian itself. A risk number that differs between Windows and Linux CI is a bug; this is where
// it would come from.

#include <cstdint>
#include <random>

namespace rl::math {

/// Deterministic seed derivation: never from time, thread id, or an address.
struct SeedPolicy {
    /// Mix (tradeId, scenarioId, stream) into a 64-bit seed with splitmix64. Distinct inputs give
    /// seeds with no visible correlation; (t, s, 0) and (t, s, 1) are independent streams (antithetic
    /// pairs use the same stream with a flag, not a second seed).
    static std::uint64_t derive(std::uint64_t tradeId, std::uint64_t scenarioId, std::uint64_t stream = 0) noexcept;
    /// splitmix64 finaliser (Steele, Lea & Flood 2014); one round.
    static std::uint64_t mix(std::uint64_t z) noexcept;
};

class Rng {
public:
    explicit Rng(std::uint64_t seed) : eng_(seed) {}

    std::uint64_t nextU64() { return eng_(); }

    /// Uniform on the open interval (0, 1): the top 53 bits of a draw, offset by half an ulp so that
    /// 0 and 1 are impossible (an exact 0 would send the inverse cdf to −∞).
    double uniform01() { return (static_cast<double>(eng_() >> 11) + 0.5) * 0x1.0p-53; }

    /// N(0,1) via the inverse cdf (AS241) — one uniform per gaussian, monotone in the uniform, which
    /// is what the Brownian bridge and Sobol need. Defined in rng.cpp to keep normal.hpp out of here.
    double gaussian();

    /// Discard n draws (for stream splitting in tests).
    void discard(std::uint64_t n) { eng_.discard(n); }

    std::mt19937_64& engine() noexcept { return eng_; }

private:
    std::mt19937_64 eng_;
};

} // namespace rl::math
