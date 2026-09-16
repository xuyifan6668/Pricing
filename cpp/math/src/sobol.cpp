#include "rl/math/sobol.hpp"

#include "rl/core/error.hpp"

#include <cstdlib>
#include <fstream>
#include <mutex>
#include <sstream>
#include <string>

#ifndef RL_DATA_DIR
#define RL_DATA_DIR "data"
#endif

namespace rl::math {

// ------------------------------------------------------------------------------------------------
// Plumbing
// ------------------------------------------------------------------------------------------------

DirectionNumbers DirectionNumbers::load(const std::filesystem::path& file) {
    std::ifstream in(file);
    if (!in) throw rl::InputError("cannot open Sobol direction-number file: " + file.string());
    DirectionNumbers dn;
    std::string line;
    std::getline(in, line); // header "d       s       a       m_i"
    RL_REQUIRE(line.find('d') != std::string::npos && line.find('m') != std::string::npos,
               "unexpected Joe-Kuo header in " + file.string());
    unsigned expected = 2;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        std::istringstream is(line);
        Dim dim;
        is >> dim.d >> dim.s >> dim.a;
        RL_REQUIRE(static_cast<bool>(is), "bad row in Sobol file at dimension " + std::to_string(expected));
        RL_REQUIRE(dim.d == expected, "Sobol file rows must be consecutive from dimension 2");
        dim.m.resize(dim.s);
        for (unsigned k = 0; k < dim.s; ++k) {
            is >> dim.m[k];
            RL_REQUIRE(static_cast<bool>(is) && (dim.m[k] % 2 == 1), "Sobol m_k must be odd integers");
        }
        dn.dims.push_back(std::move(dim));
        ++expected;
    }
    RL_REQUIRE(!dn.dims.empty(), "Sobol file has no dimensions");
    return dn;
}

const DirectionNumbers& DirectionNumbers::joeKuo() {
    static const DirectionNumbers dn = [] {
        std::filesystem::path dir = RL_DATA_DIR;
        if (const char* env = std::getenv("RL_DATA_DIR"); env && *env) dir = env;
        return load(dir / "sobol" / "new-joe-kuo-6.21201.txt");
    }();
    return dn;
}

Sobol::Sobol(unsigned dim, std::uint64_t skipCount, const DirectionNumbers& dn) {
    RL_REQUIRE(dim >= 1, "Sobol dimension must be >= 1");
    RL_REQUIRE(dim <= dn.maxDimension(), "Sobol dimension " + std::to_string(dim) + " exceeds the direction-number table (" +
                                             std::to_string(dn.maxDimension()) + ")");
    v_.assign(dim, std::vector<std::uint32_t>(kBits, 0));
    x_.assign(dim, 0);
    init(dn);
    n_ = 0;
    if (skipCount > 0) skip(skipCount);
}

void Sobol::skip(std::uint64_t n) {
    RL_REQUIRE(n_ + n < (std::uint64_t(1) << kBits), "Sobol sequence exhausted (2^32 points)");
    jumpTo(n_ + n);
}

// ------------------------------------------------------------------------------------------------
// Direction numbers and the Gray-code generator (Joe & Kuo sobol.cc, Bratley &
// Fox Algorithm 659, Antonov & Saleev 1979)
// ------------------------------------------------------------------------------------------------

void Sobol::init(const DirectionNumbers& dn) {
    const unsigned dim = dimension();
    // Dimension 0: all m_k = 1 (van der Corput in base 2).
    for (unsigned k = 0; k < kBits; ++k) v_[0][k] = 1u << (kBits - 1 - k);
    for (unsigned j = 1; j < dim; ++j) {
        const DirectionNumbers::Dim& row = dn.dims[j - 1];
        const unsigned s = row.s, a = row.a;
        std::vector<std::uint32_t>& v = v_[j];
        for (unsigned k = 0; k < kBits; ++k) {
            if (k < s) {
                v[k] = row.m[k] << (kBits - 1 - k);
            } else {
                // Recurrence from the primitive polynomial x^s + a_1 x^{s-1} + ... + a_{s-1} x + 1:
                //   v_k = v_{k-s} ^ (v_{k-s} >> s) ^ XOR_{i=1}^{s-1} [a_i v_{k-i}],  a_i = bit (s-1-i) of a.
                std::uint32_t x = v[k - s] ^ (v[k - s] >> s);
                for (unsigned i = 1; i < s; ++i)
                    if ((a >> (s - 1 - i)) & 1u) x ^= v[k - i];
                v[k] = x;
            }
        }
    }
}

std::uint64_t Sobol::next(std::span<double> out) {
    const unsigned dim = dimension();
    RL_REQUIRE(out.size() >= dim, "Sobol::next: output span too small");
    if (n_ >= (std::uint64_t(1) << kBits)) throw NumericalError("Sobol sequence exhausted (2^32 points)");
    constexpr double scale = 1.0 / 4294967296.0; // 2^-32
    for (unsigned j = 0; j < dim; ++j) out[j] = static_cast<double>(x_[j]) * scale;
    const std::uint64_t emitted = n_;
    // Gray-code advance: flip the direction number of the lowest zero bit of n_.
    unsigned c = 0;
    for (std::uint64_t m = n_; m & 1u; m >>= 1) ++c;
    if (c < kBits)
        for (unsigned j = 0; j < dim; ++j) x_[j] ^= v_[j][c];
    ++n_;
    return emitted;
}

void Sobol::jumpTo(std::uint64_t n) {
    const std::uint64_t g = n ^ (n >> 1); // Gray code of n
    for (unsigned j = 0; j < dimension(); ++j) {
        std::uint32_t x = 0;
        for (unsigned k = 0; k < kBits; ++k)
            if ((g >> k) & 1u) x ^= v_[j][k];
        x_[j] = x;
    }
    n_ = n;
}

} // namespace rl::math
