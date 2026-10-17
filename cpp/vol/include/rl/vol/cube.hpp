#pragma once
// Swaption cube queries. The cube's data and interpolation are market/swaption_cube.hpp; this
// header turns an interpolated node into a SABR smile at the caller's forward.
//
//   cubeParams(cube, e, t, F):  read ATM vol, ρ, ν at (e, t) (bilinear, total-variance ATM, flat extrapolation),
//                               then α = alphaFromAtmNormalVol(ATM, F, e, β, ρ, ν, s).
//   cubeNormalVol(..., K, F):   sabrNormalVol(cubeParams(...), F, K, e).
//
// DESIGN: the forward is an argument, not something the cube knows. The cube is market data; the forward swap
// rate of an arbitrary (expiry, tenor) is a curve calculation (linear::parRate on the underlying), which lives in
// models/ . Keeping vol independent of linear is what the §1 DAG says (math → vol, market → vol only),
// and it means a curve bump moves the smile through α's re-solve (the SABR backbone) with no cube rebuild.
//
// DESIGN: α is re-solved, never interpolated: the cube then reprices the ATM matrix exactly at any forward, which
// is the liquid, hedged part of the surface, and the interpolated (ρ, ν) shape the wings.
//
// the calibration calls cubeParams once per strip swaption and caches the result; nothing here is per path.

#include "rl/market/sabr_params.hpp"
#include "rl/market/swaption_cube.hpp"
#include "rl/math/scalar.hpp"

#include <vector>

namespace rl::vol {

/// Full SABR parameters at (expiry, tenor) in years, with α consistent with the interpolated ATM vol at `forward`.
/// Throws as SwaptionCube::interpolate and alphaFromAtmNormalVol.
template <math::Scalar T>
market::SabrParams<T> cubeParams(const market::SwaptionCube<T>& cube, double expiry, double tenor, const T& forward);

/// SABR normal vol at `strike` from the cube.
template <math::Scalar T>
T cubeNormalVol(const market::SwaptionCube<T>& cube, double expiry, double tenor, double strike, const T& forward);

/// The interpolated ATM normal vol (independent of the forward by construction).
template <math::Scalar T>
T cubeAtmNormalVol(const market::SwaptionCube<T>& cube, double expiry, double tenor);

/// Normal vols at several strikes, sharing one cubeParams call.
template <math::Scalar T>
std::vector<T> cubeSmile(const market::SwaptionCube<T>& cube, double expiry, double tenor, const T& forward, const std::vector<double>& strikes);

} // namespace rl::vol
