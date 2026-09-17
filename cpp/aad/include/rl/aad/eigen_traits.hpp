#pragma once
// Makes Eigen::Matrix<rl::aad::Real, ...> usable for products, sums and transposes — everything the
// linear-algebra-flavoured pricing code needs (Jacobian-vector products, regression design matrices).
// Decompositions (LDLT, QR) additionally need relational operators, which Real deliberately lacks
// (scalar.hpp): solvers run on double and adjoints cross them via the implicit-function theorem.

#include "rl/aad/real.hpp"

#include <Eigen/Core>

namespace Eigen {

template <>
struct NumTraits<rl::aad::Real> : NumTraits<double> {
    using Real = rl::aad::Real;
    using NonInteger = rl::aad::Real;
    using Literal = double;
    using Nested = rl::aad::Real;
    enum {
        IsComplex = 0,
        IsInteger = 0,
        IsSigned = 1,
        RequireInitialization = 1,
        ReadCost = 1,
        AddCost = 3,
        MulCost = 3,
    };
    static inline rl::aad::Real epsilon() { return NumTraits<double>::epsilon(); }
    static inline rl::aad::Real dummy_precision() { return NumTraits<double>::dummy_precision(); }
    static inline rl::aad::Real highest() { return NumTraits<double>::highest(); }
    static inline rl::aad::Real lowest() { return NumTraits<double>::lowest(); }
    static inline int digits10() { return NumTraits<double>::digits10(); }
};

} // namespace Eigen

namespace rl::aad {
using VectorR = Eigen::Matrix<Real, Eigen::Dynamic, 1>;
using MatrixR = Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>;
} // namespace rl::aad
