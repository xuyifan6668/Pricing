#include "rl/market/hash.hpp"

#include <bit>
#include <cstdio>
#include <string>

namespace rl::market {

std::uint64_t fnv1a64(const void* data, std::size_t n) noexcept {
    const auto* p = static_cast<const unsigned char*>(data);
    std::uint64_t h = Hasher::kOffsetBasis;
    for (std::size_t i = 0; i < n; ++i) {
        h ^= static_cast<std::uint64_t>(p[i]);
        h *= Hasher::kPrime;
    }
    return h;
}

Hasher& Hasher::mixBytes(const void* data, std::size_t n) noexcept {
    const auto* p = static_cast<const unsigned char*>(data);
    for (std::size_t i = 0; i < n; ++i) {
        h_ ^= static_cast<std::uint64_t>(p[i]);
        h_ *= kPrime;
    }
    return *this;
}

Hasher& Hasher::mix(std::uint64_t v) noexcept {
    // Little-endian byte order regardless of host, so the hash is platform-independent.
    unsigned char bytes[8];
    for (int i = 0; i < 8; ++i) bytes[i] = static_cast<unsigned char>((v >> (8 * i)) & 0xFFu);
    return mixBytes(bytes, 8);
}

Hasher& Hasher::mix(double v) noexcept {
    if (v == 0.0) v = 0.0; // canonicalise −0.0
    return mix(std::bit_cast<std::uint64_t>(v));
}

Hasher& Hasher::mix(std::string_view s) noexcept {
    mix(static_cast<std::uint64_t>(s.size()));
    return mixBytes(s.data(), s.size());
}

Hasher& Hasher::mix(const std::vector<double>& v) noexcept {
    mix(static_cast<std::uint64_t>(v.size()));
    for (double x : v) mix(x);
    return *this;
}

std::string toHex(std::uint64_t h) {
    char buf[2 + 16 + 1];
    std::snprintf(buf, sizeof buf, "0x%016llx", static_cast<unsigned long long>(h));
    return std::string(buf);
}

} // namespace rl::market
