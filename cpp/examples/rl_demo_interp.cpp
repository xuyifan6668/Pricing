// rl_demo_interp NPOINTS t0:y0 t1:y1 ...   (y = zero rate; nodes at increasing t, first t must be 0)
//   e.g. rl_demo_interp 400 0:0 0.5:0.02 1:0.024 2:0.03 5:0.039 10:0.045 30:0.052
// Prints CSV: t, linear_y, linear_f, mc_y, mc_f, cubic_y, cubic_f  where y = -log P(0,t) = r·t and
// f = dy/dt is the instantaneous forward under each interpolator.
#include "rl/core/error.hpp"
#include "rl/math/interpolation.hpp"

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

int main(int argc, char** argv) {
    using namespace rl::math;
    try {
        if (argc < 4) {
            std::cerr << "usage: rl_demo_interp NPOINTS t0:r0 t1:r1 ...\n";
            return 2;
        }
        const int n = std::atoi(argv[1]);
        std::vector<double> ts, ys;
        for (int i = 2; i < argc; ++i) {
            std::string a = argv[i];
            auto colon = a.find(':');
            RL_REQUIRE(colon != std::string::npos, "node must be t:r");
            double t = std::stod(a.substr(0, colon));
            double r = std::stod(a.substr(colon + 1));
            ts.push_back(t);
            ys.push_back(r * t);
        }
        Linear<double> lin(ts, ys);
        MonotoneConvex<double> mc(ts, ys);
        NaturalCubic<double> cub(ts, ys, Extrapolation::Linear);
        const double tmax = ts.back() * 1.1;
        std::printf("t,linear_y,linear_f,mc_y,mc_f,cubic_y,cubic_f\n");
        for (int i = 0; i <= n; ++i) {
            double t = tmax * i / n;
            std::printf("%.6f,%.12f,%.12f,%.12f,%.12f,%.12f,%.12f\n", t, lin(t), lin.derivative(t), mc(t), mc.forward(t), cub(t), cub.derivative(t));
        }
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << "\n";
        return 1;
    }
}
