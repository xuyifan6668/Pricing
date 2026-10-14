#pragma once
// Normal ↔ shifted-lognormal vol conversion.
//
//   normalToShiftedLognormalVol, shiftedLognormalToNormalVol — the exact, price-matching
//   conversion and the shift handling (reasoning in src/conversion.cpp).
//   atmNormalFromLognormalApprox, atmLognormalFromNormalApprox — the ATM series, for the
//   notebook and as an initial guess.
//
// Why a desk cares: every EUR/CHF/JPY screen quotes shifted-lognormal vols with a *published* shift and every
// USD screen quotes normal vols; models calibrate on one of them. A conversion that is not exact at the 1e-12
// level turns into a calibration residual that nobody can explain.
//
// double-only: conversion runs once per quote at cube build time, never inside a pricing, so nothing is
// differentiated through it (the cube stores normal vols).

namespace rl::vol {

/// σ_LN with shift s such that the shifted-Black price at σ_LN equals the Bachelier price at σ_N, for the same
/// forward, strike and expiry. Independent of call/put by parity. See src/conversion.cpp for the spec.
double normalToShiftedLognormalVol(double forward, double strike, double expiry, double normalVol, double shift);

/// σ_N such that the Bachelier price at σ_N equals the shifted-Black price at σ_LN with shift s.
double shiftedLognormalToNormalVol(double forward, double strike, double expiry, double lognormalVol, double shift);

/// ATM series from equating the ATM prices: σ_N √T/√(2π) = (F+s)(2Φ(σ_LN√T/2) − 1)
///   ⇒ σ_N = σ_LN (F + s) (1 − σ_LN² T/24 + σ_LN⁴ T²/640 − …); this returns the first two terms.
/// Relative error ≈ σ_LN⁴T²/640 (4e-4 at σ√T = 1). Throws InputError if F + s ≤ 0.
double atmNormalFromLognormalApprox(double forward, double expiry, double lognormalVol, double shift);

/// Inverse of the above to the same order: σ_LN = σ_N/(F + s) · (1 + σ_N² T / (24 (F + s)²)).
double atmLognormalFromNormalApprox(double forward, double expiry, double normalVol, double shift);

} // namespace rl::vol
