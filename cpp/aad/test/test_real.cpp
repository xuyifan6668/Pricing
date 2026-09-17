#include "rl/aad/real.hpp"

#include "rl/math/scalar.hpp"

#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <utility>
#include <vector>
#include <sstream>

using rl::aad::Real;
using rl::aad::Tape;

static_assert(rl::math::Scalar<Real>, "Real must satisfy the Scalar concept");

TEST_CASE("Constants are passive and record nothing", "[aad][real]") {
    Tape& t = Tape::current();
    t.reset();
    Real a = 2.0;
    Real b(3.0);
    Real c = a * b + exp(a) / b - sqrt(b);
    CHECK_FALSE(c.isActive());
    CHECK(c.value() == 2.0 * 3.0 + std::exp(2.0) / 3.0 - std::sqrt(3.0));
    CHECK(t.empty());
}

TEST_CASE("Active leaves and every operator record one node per operation", "[aad][real]") {
    Tape& t = Tape::current();
    t.reset();
    Real x = Real::input(1.5);
    Real y = Real::input(-0.25);
    REQUIRE(t.size() == 2);
    Real z = x + y;   // 3
    z = z - y;        // 4
    z = z * x;        // 5
    z = z / y;        // 6
    z = z + 1.0;      // 7
    z = 2.0 - z;      // 8
    z = z * 3.0;      // 9
    z = 4.0 / z;      // 10
    z = -z;           // 11
    z = +z;           // no node
    z = exp(z);       // 12
    z = log(fabs(z)); // 13, 14
    z = sqrt(fabs(z));// 15, 16
    z = erf(z) + erfc(z); // 17, 18, 19
    z = pow(z, 2.0);  // 20
    z = pow(z, x);    // 21
    z = pow(2.0, z);  // 22
    z += x;           // 23
    z -= y;           // 24
    z *= x;           // 25
    z /= y;           // 26
    z += 1.0;         // 27
    z -= 1.0;         // 28
    z *= 2.0;         // 29
    z /= 2.0;         // 30
    CHECK(t.size() == 30);
    CHECK(z.isActive());
    CHECK(z.node() == 29);
    // Two-parent nodes have arity 2 when both parents are active.
    CHECK(t.arity(2) == 2);
    CHECK(t.arity(6) == 1); // z + 1.0
    t.reset();
}

TEST_CASE("Active arithmetic reproduces double bit for bit", "[aad][real]") {
    Tape& t = Tape::current();
    t.reset();
    const double xv = 0.7, yv = 1.9;
    Real x = Real::input(xv), y = Real::input(yv);
    auto f = [](auto a, auto b) {
        using std::exp;
        using std::log;
        using std::sqrt;
        using std::pow;
        using std::erf;
        return (a * b - a / b + 3.0 * a) * exp(-b) + log(a + b) - sqrt(a * a + b * b) + pow(a, b) - erf(a - b) + pow(a, 1.5) + pow(2.0, b);
    };
    Real r = f(x, y);
    double d = f(xv, yv);
    CHECK(r.value() == d);
    CHECK(r.isActive());
    t.reset();
}

TEST_CASE("Mixed constant/active operations keep only the active parent", "[aad][real]") {
    Tape& t = Tape::current();
    t.reset();
    Real x = Real::input(2.0);
    Real k = 5.0;
    Real z = x * k;
    CHECK(t.arity(z.node()) == 1);
    CHECK(t.parent(z.node(), 0) == x.node());
    Real w = k + k; // constant + constant: no node
    CHECK_FALSE(w.isActive());
    CHECK(t.size() == 2);
    t.reset();
}

TEST_CASE("max/min/positivePart return one of the inputs (identity on the tape)", "[aad][real]") {
    Tape& t = Tape::current();
    t.reset();
    Real x = Real::input(2.0), y = Real::input(3.0);
    const Real& m = rl::math::max(x, y);
    CHECK(&m == &y);
    CHECK(&rl::math::min(x, y) == &x);
    Real neg = Real::input(-1.0);
    Real p = rl::math::positivePart(neg);
    CHECK_FALSE(p.isActive());
    CHECK(p.value() == 0.0);
    Real q = rl::math::positivePart(x);
    CHECK(q.node() == x.node());
    CHECK(rl::math::sign(neg) == -1.0);
    CHECK(rl::math::isFinite(x));
    CHECK(t.size() == 3);
    t.reset();
}

namespace {
// One function touching every operator and elementary function Real overloads.
template <class T>
T everything(const T& x, const T& y) {
    using std::erf; using std::erfc; using std::exp; using std::fabs; using std::log; using std::pow; using std::sqrt;
    T a = x * y + x / y - (x - y) + 2.0 * x - y / 3.0 + 1.5 / x + (4.0 - y);
    T b = exp(-x) * log(y + 2.0) + sqrt(x * x + 1.0) + fabs(x - 2.0 * y);
    T c = erf(x / y) - erfc(y) + pow(x, y) + pow(x, 1.7) + pow(1.3, y);
    T d = -a;
    d += b;
    d *= c;
    d -= x;
    d /= y;
    return d;
}
} // namespace

TEST_CASE("Adjoints of every operator match central differences", "[aad][real][adjoint]") {
    Tape& t = Tape::current();
    for (auto [x0, y0] : {std::pair{0.7, 1.9}, std::pair{1.3, 0.4}, std::pair{2.2, 2.9}}) {
        t.reset();
        Real x = Real::input(x0), y = Real::input(y0);
        Real f = everything(x, y);
        CHECK(f.value() == everything(x0, y0));
        t.seedAdjoint(f.node(), 1.0);
        t.propagate();
        const double h = 1e-6;
        const double fx = (everything(x0 + h, y0) - everything(x0 - h, y0)) / (2.0 * h);
        const double fy = (everything(x0, y0 + h) - everything(x0, y0 - h)) / (2.0 * h);
        CHECK(std::fabs(t.adjoint(x.node()) - fx) < 1e-7 * (1.0 + std::fabs(fx)));
        CHECK(std::fabs(t.adjoint(y.node()) - fy) < 1e-7 * (1.0 + std::fabs(fy)));
    }
    t.reset();
}

TEST_CASE("One sweep gives the whole gradient: f = sum of x_i^2 over 1000 inputs", "[aad][real][adjoint]") {
    Tape& t = Tape::current();
    t.reset();
    std::vector<Real> xs;
    for (int i = 0; i < 1000; ++i) xs.push_back(Real::input(0.001 * i));
    Real f = 0.0;
    for (const auto& x : xs) f += x * x;
    t.seedAdjoint(f.node(), 1.0);
    t.propagate();
    for (int i = 0; i < 1000; ++i) CHECK(t.adjoint(xs[static_cast<std::size_t>(i)].node()) == 2.0 * 0.001 * i);
    t.reset();
}

TEST_CASE("Streaming", "[aad][real]") {
    Tape& t = Tape::current();
    t.reset();
    std::ostringstream os;
    os << Real(1.5) << " " << Real::input(2.5);
    CHECK(os.str() == "1.5 2.5[#0]");
    t.reset();
}
