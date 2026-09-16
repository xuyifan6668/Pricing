#include "rl/math/rng.hpp"

#include "rl/math/normal.hpp"

namespace rl::math {

std::uint64_t SeedPolicy::mix(std::uint64_t z) noexcept {
    z += 0x9E3779B97F4A7C15ull;
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
    return z ^ (z >> 31);
}

std::uint64_t SeedPolicy::derive(std::uint64_t tradeId, std::uint64_t scenarioId, std::uint64_t stream) noexcept {
    // Three rounds so that each input is fully diffused before the next is folded in.
    std::uint64_t s = mix(tradeId);
    s = mix(s ^ (scenarioId + 0x632BE59BD9B4E019ull));
    s = mix(s ^ (stream + 0x9E6C63D0676A9A99ull));
    return s;
}

double Rng::gaussian() { return normalInverseCdf<double>(uniform01()); }

} // namespace rl::math
