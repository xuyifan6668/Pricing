#pragma once
// definitions of the cube queries (cube.hpp).
#include "rl/vol/cube.hpp"
#include "rl/vol/sabr.hpp"

namespace rl::vol {

template <math::Scalar T>
market::SabrParams<T> cubeParams(const market::SwaptionCube<T>& cube, double expiry, double tenor, const T& forward) {
    const market::InterpolatedNode<T> in = cube.interpolate(expiry, tenor);
    market::SabrParams<T> p;
    p.rho = in.rho;
    p.nu = in.nu;
    p.beta = in.beta;
    p.shift = in.shift;
    p.alpha = alphaFromAtmNormalVol(in.atmNormalVol, forward, expiry, in.beta, in.rho, in.nu, in.shift);
    return p;
}

template <math::Scalar T>
T cubeNormalVol(const market::SwaptionCube<T>& cube, double expiry, double tenor, double strike, const T& forward) {
    return sabrNormalVol(cubeParams(cube, expiry, tenor, forward), forward, strike, expiry);
}

template <math::Scalar T>
T cubeAtmNormalVol(const market::SwaptionCube<T>& cube, double expiry, double tenor) { return cube.interpolate(expiry, tenor).atmNormalVol; }

template <math::Scalar T>
std::vector<T> cubeSmile(const market::SwaptionCube<T>& cube, double expiry, double tenor, const T& forward, const std::vector<double>& strikes) {
    const market::SabrParams<T> p = cubeParams(cube, expiry, tenor, forward);
    std::vector<T> out;
    out.reserve(strikes.size());
    for (double k : strikes) out.push_back(sabrNormalVol(p, forward, k, expiry));
    return out;
}

} // namespace rl::vol

#define RL_VOL_INSTANTIATE_CUBE(T)                                                                                          \
    namespace rl::vol {                                                                                                     \
    template market::SabrParams<T> cubeParams<T>(const market::SwaptionCube<T>&, double, double, const T&);                 \
    template T cubeNormalVol<T>(const market::SwaptionCube<T>&, double, double, double, const T&);                           \
    template T cubeAtmNormalVol<T>(const market::SwaptionCube<T>&, double, double);                                         \
    template std::vector<T> cubeSmile<T>(const market::SwaptionCube<T>&, double, double, const T&, const std::vector<double>&); \
    }
