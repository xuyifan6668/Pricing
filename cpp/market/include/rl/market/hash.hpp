#pragma once
// Content hashing for the market snapshot: the hash goes on every pricing result and every
// report line, so it must be a pure function of the market *content*, identical on MSVC and GCC, and
// sensitive to a one-ulp change of any input.
//
// FNV-1a, 64-bit (Fowler–Noll–Vo; public domain; see Noll's reference page). Not cryptographic and not
// meant to be: it identifies a snapshot, it does not authenticate one. Doubles are mixed by their IEEE-754
// bit pattern (std::bit_cast), never by formatting; +0.0 and −0.0 are canonicalised to +0.0 so that a
// bump of exactly zero is not a "different market". NaN is mixed as its bit pattern and will make two
// otherwise-equal snapshots differ, which is the loud behaviour we want.
//
// Strings are length-prefixed so that ("ab","c") and ("a","bc") hash differently.

#include "rl/core/date.hpp"

#include <cstddef>
#include <cstdint>
#include <string_view>
#include <vector>

namespace rl::market {

class Hasher {
public:
    static constexpr std::uint64_t kOffsetBasis = 14695981039346656037ull;
    static constexpr std::uint64_t kPrime = 1099511628211ull;

    constexpr Hasher() noexcept = default;

    Hasher& mixBytes(const void* data, std::size_t n) noexcept;
    Hasher& mix(std::uint64_t v) noexcept;
    Hasher& mix(std::int64_t v) noexcept { return mix(static_cast<std::uint64_t>(v)); }
    Hasher& mix(std::uint32_t v) noexcept { return mix(static_cast<std::uint64_t>(v)); }
    Hasher& mix(int v) noexcept { return mix(static_cast<std::int64_t>(v)); }
    Hasher& mix(bool v) noexcept { return mix(static_cast<std::uint64_t>(v ? 1 : 0)); }
    Hasher& mix(double v) noexcept;
    Hasher& mix(std::string_view s) noexcept;
    Hasher& mix(const char* s) noexcept { return mix(std::string_view(s)); }
    Hasher& mix(core::Date d) noexcept { return mix(static_cast<std::int64_t>(d.serial())); }
    Hasher& mix(const std::vector<double>& v) noexcept;

    constexpr std::uint64_t value() const noexcept { return h_; }

private:
    std::uint64_t h_ = kOffsetBasis;
};

/// FNV-1a of a raw byte range (what Hasher::mixBytes does from the offset basis).
std::uint64_t fnv1a64(const void* data, std::size_t n) noexcept;

/// "0x" + 16 lowercase hex digits; the form printed in reports.
std::string toHex(std::uint64_t h);

} // namespace rl::market
