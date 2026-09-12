// standard one-dimensional root finders, transcribed from Numerical Recipes
// 3e (§9.1 zbrac, §9.3 zbrent, §9.4 rtsafe) with the error policy of rootfind.hpp.
#include "rl/math/rootfind.hpp"

#include "rl/core/error.hpp"

#include <cmath>
#include <limits>

namespace rl::math {

RootResult brent(const std::function<double(double)>& f, double a, double b, RootOptions opt) {
    RL_REQUIRE(opt.maxIter > 0 && opt.xTol >= 0.0, "brent: bad options");
    double fa = f(a), fb = f(b);
    if (fa == 0.0) return RootResult{a, 0.0, 0, true};
    if (fb == 0.0) return RootResult{b, 0.0, 0, true};
    if ((fa > 0.0) == (fb > 0.0)) throw InputError("brent: root is not bracketed");
    double c = a, fc = fa, d = b - a, e = d;
    constexpr double eps = std::numeric_limits<double>::epsilon();
    for (int iter = 1; iter <= opt.maxIter; ++iter) {
        if ((fb > 0.0) == (fc > 0.0)) { // c must bracket with b
            c = a;
            fc = fa;
            d = e = b - a;
        }
        if (std::fabs(fc) < std::fabs(fb)) { // b is the best estimate
            a = b; b = c; c = a;
            fa = fb; fb = fc; fc = fa;
        }
        const double tol1 = 2.0 * eps * std::fabs(b) + 0.5 * opt.xTol;
        const double m = 0.5 * (c - b);
        if (std::fabs(m) <= tol1 || fb == 0.0 || (opt.fTol > 0.0 && std::fabs(fb) <= opt.fTol)) return RootResult{b, fb, iter, true};
        if (std::fabs(e) >= tol1 && std::fabs(fa) > std::fabs(fb)) {
            const double s = fb / fa;
            double p, q;
            if (a == c) { // secant
                p = 2.0 * m * s;
                q = 1.0 - s;
            } else { // inverse quadratic interpolation
                const double qa = fa / fc, r = fb / fc;
                p = s * (2.0 * m * qa * (qa - r) - (b - a) * (r - 1.0));
                q = (qa - 1.0) * (r - 1.0) * (s - 1.0);
            }
            if (p > 0.0) q = -q;
            p = std::fabs(p);
            const double min1 = 3.0 * m * q - std::fabs(tol1 * q), min2 = std::fabs(e * q);
            if (2.0 * p < std::min(min1, min2)) { // accept interpolation
                e = d;
                d = p / q;
            } else { // bisection
                d = m;
                e = m;
            }
        } else {
            d = m;
            e = m;
        }
        a = b;
        fa = fb;
        b += std::fabs(d) > tol1 ? d : (m > 0.0 ? tol1 : -tol1);
        fb = f(b);
    }
    throw NumericalError("brent: no convergence in " + std::to_string(opt.maxIter) + " iterations");
}

RootResult newtonSafeguarded(const std::function<std::pair<double, double>(double)>& fdf, double x0, double lo, double hi, RootOptions opt) {
    RL_REQUIRE(lo < hi, "newtonSafeguarded: need lo < hi");
    RL_REQUIRE(x0 >= lo && x0 <= hi, "newtonSafeguarded: x0 outside [lo, hi]");
    auto [flo, dlo] = fdf(lo);
    auto [fhi, dhi] = fdf(hi);
    if (flo == 0.0) return RootResult{lo, 0.0, 0, true};
    if (fhi == 0.0) return RootResult{hi, 0.0, 0, true};
    RL_REQUIRE((flo > 0.0) != (fhi > 0.0), "newtonSafeguarded: root is not bracketed");
    // Orient so that f(xl) < 0 < f(xh).
    double xl = lo, xh = hi;
    if (flo > 0.0) std::swap(xl, xh);
    double x = x0;
    auto [fx, dfx] = fdf(x);
    double dxold = std::fabs(hi - lo), dx = dxold;
    for (int iter = 1; iter <= opt.maxIter; ++iter) {
        const bool newtonOut = ((x - xh) * dfx - fx) * ((x - xl) * dfx - fx) > 0.0; // Newton step would leave the bracket
        const bool tooSlow = std::fabs(2.0 * fx) > std::fabs(dxold * dfx);           // not halving fast enough
        if (newtonOut || tooSlow || dfx == 0.0) {
            dxold = dx;
            dx = 0.5 * (xh - xl);
            x = xl + dx;
        } else {
            dxold = dx;
            dx = fx / dfx;
            x -= dx;
        }
        std::tie(fx, dfx) = fdf(x);
        if (std::fabs(dx) < opt.xTol || fx == 0.0 || (opt.fTol > 0.0 && std::fabs(fx) <= opt.fTol)) return RootResult{x, fx, iter, true};
        if (fx < 0.0) xl = x;
        else xh = x;
    }
    throw NumericalError("newtonSafeguarded: no convergence in " + std::to_string(opt.maxIter) + " iterations");
}

bool expandBracket(const std::function<double(double)>& f, double& a, double& b, int maxExpansions) {
    RL_REQUIRE(a != b, "expandBracket: a == b");
    constexpr double factor = 1.6;
    double fa = f(a), fb = f(b);
    for (int i = 0; i < maxExpansions; ++i) {
        if ((fa > 0.0) != (fb > 0.0) || fa == 0.0 || fb == 0.0) return true;
        if (std::fabs(fa) < std::fabs(fb)) fa = f(a += factor * (a - b)); // move the endpoint nearer the root
        else fb = f(b += factor * (b - a));
    }
    return (fa > 0.0) != (fb > 0.0);
}

} // namespace rl::math
