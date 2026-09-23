# rl::market

Quotes, fixings, the content hash, and the curve / curve-set / snapshot templates.

## Contract

- `QuoteSet` is a validated, immutable list of `Quote`s with an as-of date; JSON in `data/quotes/`.
  `withValues`/`bumped` return copies. A quote's `curve` field says which curve the instrument pins
  (default: its index; a basis swap may be assigned to the other curve — that is what couples the curves).
- `FixingStore` holds realised fixings as passive `double`, keyed `(index name, fixing date)`;
  the fixing *for* date `d` is known from `d + 1` business day onwards. JSON in `data/fixings/`, percent.
- `Hasher` is FNV-1a 64 over IEEE bit patterns, length-prefixed strings, little-endian integers:
  identical on MSVC and GCC; `−0.0` canonicalised to `+0.0`.
- `CurveInterpolation` names the interpolator applied to `y = −log P(0,t)`: `LinearLogDf` (piecewise-flat
  forwards) or `MonotoneConvex` (Hagan–West).
- `YieldCurve<T>`: pillars in ACT/365F years, ordinates `y = −log P`, implicit `t = 0` anchor,
  interpolator variant selected by `CurveInterpolation`, flat-forward extrapolation both ways, and an
  additive turn-of-year jump layer (`turns()`, `withTurnJumps`). `df(0) == T(1)` exactly.
- `CurveSet<T>`: named curves plus a discount role and an index→curve projection map; the SOFR curve
  is stored once. `withCurve` replaces one curve and keeps the roles.
- `MarketSnapshot<T>`: no non-const member, private constructor, built only by `SnapshotBuilder<T>`
  into `shared_ptr<const>`; hash computed once over as-of, curves (via `value()`), fixings, quotes, and the swaption
  cube when present. `swaptionCube()` is `nullptr` for a curves-only snapshot.
- `SabrParams<T>`: α, ρ, ν active; β and shift passive policy; `validateSabrParams` checks the domain.
- `SwaptionVolQuotes`: the vol screen — ATM and smile offsets per expiry × tenor in normal vols, with the node forward
  as a `// SIMPLIFICATION:`; a complete grid is required; JSON in `data/vol/` in human units, decimals in memory.
- `SwaptionCube<T>`: expiry × tenor grid of `CubeNode<T>` (forward, market ATM vol, fitted SABR). `interpolate` is
  bilinear with ATM as total variance and flat extrapolation; α is not interpolated (`vol::cubeParams` re-solves it). JSON for `double`.

## Invariants (tested)

- `Quote::validate` rejects missing dates/tenors, wrong units (a futures price as a rate), |rate| ≥ 100%.
- `QuoteSet` rejects duplicate ids; JSON round-trips; errors carry the quote id.
- One ulp on any double changes the hash; equal content gives equal hashes; `Real` and `double`
  snapshots of the same market hash equally.
- `YieldCurve`: node reproduction, `df(0) == 1`, positivity, `Real == double` bit for bit,
  `withNodeValues` leaves the original untouched and validates its input.
- Cube: nodes reproduced, total variance linear in expiry, weights sum to one, flat beyond both axes, JSON round trip
  in any node order and loud on holes and duplicates, one ulp changes the hash, `Real` hashes as `double`.
- Vol quotes: holes, duplicates, zero or unsorted offsets and unit errors are rejected with the node id; the shipped
  file is a complete 9 × 9 grid with positive strikes.
- Snapshot: a cube changes the hash, a bumped cube changes it again, and a cube with another as-of is rejected.

## Threading

Everything is immutable after construction; `FixingStore` and `QuoteSet` are value types mutated only
while a snapshot is assembled. Share `shared_ptr<const MarketSnapshot<T>>` freely across the risk driver.
