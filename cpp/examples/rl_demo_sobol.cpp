// rl_demo_sobol points N DIM        -> CSV of N points (Sobol, skipping the origin) and N MT points
// rl_demo_sobol integrate KMAX DIM  -> CSV: N, sobolErr, mtErr for N = 2^4 .. 2^KMAX on exp(-|x|_1)
#include "rl/core/error.hpp"
#include "rl/math/rng.hpp"
#include "rl/math/sobol.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

int main(int argc, char** argv) {
    using namespace rl::math;
    try {
        if (argc != 4) {
            std::cerr << "usage: rl_demo_sobol points N DIM | rl_demo_sobol integrate KMAX DIM\n";
            return 2;
        }
        const std::string mode = argv[1];
        const unsigned dim = static_cast<unsigned>(std::atoi(argv[3]));
        if (mode == "points") {
            const std::size_t n = static_cast<std::size_t>(std::atoll(argv[2]));
            Sobol s(dim, 1);
            Rng rng(SeedPolicy::derive(1, 1));
            std::vector<double> p(dim);
            std::printf("generator,index");
            for (unsigned j = 0; j < dim; ++j) std::printf(",x%u", j + 1);
            std::printf("\n");
            for (std::size_t i = 0; i < n; ++i) {
                s.next(p);
                std::printf("sobol,%zu", i);
                for (double v : p) std::printf(",%.10f", v);
                std::printf("\n");
            }
            for (std::size_t i = 0; i < n; ++i) {
                std::printf("mt,%zu", i);
                for (unsigned j = 0; j < dim; ++j) std::printf(",%.10f", rng.uniform01());
                std::printf("\n");
            }
        } else if (mode == "integrate") {
            const int kmax = std::atoi(argv[2]);
            const double exact = std::pow(1.0 - std::exp(-1.0), dim);
            std::printf("N,sobolErr,mtErr\n");
            for (int k = 4; k <= kmax; ++k) {
                const std::size_t n = std::size_t(1) << k;
                Sobol s(dim, 1);
                std::vector<double> p(dim);
                double acc = 0.0;
                for (std::size_t i = 0; i < n; ++i) {
                    s.next(p);
                    double e = 0.0;
                    for (double v : p) e += v;
                    acc += std::exp(-e);
                }
                const double sobolErr = std::fabs(acc / static_cast<double>(n) - exact);
                double mtErr = 0.0;
                const int seeds = 16;
                for (int sd = 0; sd < seeds; ++sd) {
                    Rng rng(SeedPolicy::derive(42, static_cast<std::uint64_t>(sd)));
                    double a2 = 0.0;
                    for (std::size_t i = 0; i < n; ++i) {
                        double e = 0.0;
                        for (unsigned j = 0; j < dim; ++j) e += rng.uniform01();
                        a2 += std::exp(-e);
                    }
                    mtErr += std::fabs(a2 / static_cast<double>(n) - exact);
                }
                std::printf("%zu,%.10e,%.10e\n", n, sobolErr, mtErr / seeds);
            }
        } else {
            std::cerr << "unknown mode " << mode << "\n";
            return 2;
        }
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << "\n";
        return 1;
    }
}
