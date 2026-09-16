#pragma once
// Sobol low-discrepancy sequence with Joe–Kuo direction numbers (Joe & Kuo 2008, "Constructing
// Sobol sequences with better two-dimensional projections", SIAM J. Sci. Comput. 30, 2635–2654;
// file new-joe-kuo-6.21201 from https://web.maths.unsw.edu.au/~fkuo/sobol/, BSD licence).
//
// Conventions (be consistent with the Brownian bridge!):
//   * 32-bit integers; point n has coordinates xⱼ(n) = vⱼ(n) / 2³², so the first point (n = 0) is the
//     origin. That matches the Joe–Kuo reference implementation exactly, which is what the test
//     checks. Callers feeding the inverse normal cdf MUST skip the origin: construct with skip ≥ 1.
//     the path generator skips exactly one point and says so in its README.
//   * Gray-code ordering: point n is obtained from point n−1 by XOR with one direction number, so
//     generation is O(dim) per point. The points are the same *set* as the natural ordering over any
//     block of 2^m consecutive points starting at a multiple of 2^m.
//   * Dimension j (1-based in the file) uses the file's row for d = j; dimension 1 has all mₖ = 1.

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <span>
#include <vector>

namespace rl::math {

struct DirectionNumbers {
    struct Dim {
        unsigned d = 0;              ///< dimension index in the file (2, 3, ...)
        unsigned s = 0;              ///< degree of the primitive polynomial
        unsigned a = 0;              ///< polynomial coefficients (bits), per Joe–Kuo
        std::vector<unsigned> m;     ///< initial direction integers m₁..mₛ (odd)
    };
    std::vector<Dim> dims;           ///< dims[0] is file dimension 2; dimension 1 is implicit.

    std::size_t maxDimension() const noexcept { return dims.size() + 1; }
    /// Parse a Joe–Kuo file ("d s a m_i" header then one row per dimension). Throws InputError.
    static DirectionNumbers load(const std::filesystem::path& file);
    /// The bundled data/sobol/new-joe-kuo-6.21201.txt (truncated to 1111 dimensions), loaded once.
    static const DirectionNumbers& joeKuo();
};

class Sobol {
public:
    static constexpr unsigned kBits = 32;

    /// `dim` ≥ 1 dimensions; `skip` points discarded before the first call to next().
    explicit Sobol(unsigned dim, std::uint64_t skip = 0, const DirectionNumbers& dn = DirectionNumbers::joeKuo());

    unsigned dimension() const noexcept { return static_cast<unsigned>(v_.size()); }
    /// Index of the point that the next call to next() returns.
    std::uint64_t index() const noexcept { return n_; }

    /// Write the next point into `out` (size ≥ dimension()); returns the point's index.
    std::uint64_t next(std::span<double> out);
    std::vector<double> next() {
        std::vector<double> p(dimension());
        next(p);
        return p;
    }
    /// Discard n points (uses the direct formula, O(dim · 32), not n calls to next()).
    void skip(std::uint64_t n);

    /// Raw integer coordinates of the current point (tests).
    const std::vector<std::uint32_t>& state() const noexcept { return x_; }

private:
    void init(const DirectionNumbers& dn);         // the maths: v_[j][k] from (s, a, m)
    void jumpTo(std::uint64_t n);                  // the maths: set x_ to point n directly

    std::vector<std::vector<std::uint32_t>> v_;    // v_[j][k], j < dim, k < 32: direction numbers scaled to 32 bits
    std::vector<std::uint32_t> x_;                 // current point (integers), for the Gray-code update
    std::uint64_t n_ = 0;                          // index of the point x_ currently holds (returned by next())
};

} // namespace rl::math
