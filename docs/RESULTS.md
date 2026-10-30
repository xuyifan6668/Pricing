# Results

Everything here is produced by `rl_analysis` (C++) and rendered by `tools/make_results.py`; nothing is typed
in by hand. Market: USD, as of 2026-09-04, from the files in `data/` (the curve quotes and the vol screen are
synthetic but realistic). Prices are in basis points of notional unless stated.

Timings: Intel(R) Xeon(R) Processor @ 2.10GHz, 4 logical cores, Linux x86_64, single-threaded, GCC -O2.

To regenerate:

```bash
cmake -S . -B build/dev -G Ninja -DCMAKE_BUILD_TYPE=Release && cmake --build build/dev
```

```bash
uv run python tools/make_results.py
```

## 1. Curves: USD SOFR and Fed Funds, bootstrapped together

**Product view.** 37 instruments define two curves: SOFR overnight-index swaps and 3-month SOFR futures for
the SOFR curve, Fed Funds–SOFR basis swaps and Fed Funds swaps for the Fed Funds curve. The 25- and 30-year
SOFR points are quoted only through basis swaps, so the SOFR long end depends on the Fed Funds curve while
Fed Funds cashflows are discounted on SOFR. Neither curve can be built before the other.

| Interpolation | Method | Newton iterations | Worst residual (bp) | Worst instrument | Time (ms) |
| --- | --- | ---: | ---: | --- | ---: |
| linear_logdf | global Newton (adjoint Jacobian) | 4 | 4.7e-11 | SOFR-OIS-1W | 4.7 |
| linear_logdf | global Newton (finite-difference Jacobian) | 4 | 6.7e-11 | SOFR-OIS-1W | 39.6 |
| linear_logdf | sequential (one pass curve by curve) | 5 | 5.0e+00 | FF-SOFR-25Y | 5.4 |
| monotone_convex | global Newton (adjoint Jacobian) | 4 | 6.4e-11 | FF-SOFR-30Y | 11.2 |
| monotone_convex | global Newton (finite-difference Jacobian) | 4 | 4.7e-11 | SOFR-OIS-1W | 57.5 |
| monotone_convex | sequential (one pass curve by curve) | 8 | 5.0e+00 | FF-SOFR-25Y | 20.2 |

- The global solve reprices every instrument to about 1e-10 bp in four iterations under both interpolations
  (worst residuals: linear_logdf 4.7e-11 bp, monotone_convex 6.4e-11 bp).
- One pass of the curve-by-curve bootstrap leaves the long basis swaps about 5 bp off, because they were
  solved against a Fed Funds curve that did not exist yet.
- The Jacobian is taken by an adjoint sweep through the same pricing code the solver runs. It gives the same
  answer as finite differences about ten times faster, and stays exact under monotone convex interpolation,
  where a hand-derived Jacobian is easy to get subtly wrong.

![curve_forwards.png](img/curve_forwards.png)

*Interpolation is a modelling choice, not a detail.* Linear interpolation of log discount factors gives
piecewise-constant forwards with a jump at every pillar; monotone convex gives continuous forwards. Both
reprice the same quotes exactly. They differ in the forward rate between pillars, which is what a
forward-starting swap or a futures strip is priced on.

## 2. Swap risk: adjoint differentiation against bumping

**Product view.** A swap desk hedges with the instruments the curve is built from, so risk is reported per
par quote: "how much does the book move if the 10-year swap rate rises 1 bp and the curve is rebuilt".

A book of 200 swaps (1–30 years, both directions), delta to all 37 curve quotes:

| Method | Time (ms) | Relative to one valuation |
| --- | ---: | ---: |
| one valuation of the book | 0.38 | 1× |
| adjoint + implicit-function theorem | 2.23 | 6× |
| bump and re-bootstrap each quote | 184.48 | 485× |

- Adjoint: one reverse sweep gives the sensitivity to every curve node; one transposed linear solve with the
  bootstrap Jacobian (implicit-function theorem) converts node risk to par-quote risk. **83× faster**
  than bumping here, and the cost does not grow with the number of quotes.
- The two agree to 1e-05 of the largest delta when the bump is small (0.01 bp). A 1 bp bump differs
  by about 1e-3 because it also contains the gamma.

Largest deltas (currency units per +1 bp, notionals 1–5 million):

| Quote | Adjoint | Bump 1 bp | Bump 0.01 bp |
| --- | ---: | ---: | ---: |
| FF-SOFR-30Y | 101,850.4 | 101,870.0 | 101,850.6 |
| FF-OIS-30Y | 101,646.9 | 101,700.3 | 101,647.5 |
| SOFR-OIS-15Y | 54,407.3 | 54,432.0 | 54,407.5 |
| FF-OIS-25Y | 40,055.8 | 40,157.6 | 40,056.9 |
| FF-SOFR-25Y | 39,754.5 | 39,789.3 | 39,754.9 |
| SOFR-OIS-12Y | 27,846.4 | 27,846.8 | 27,846.4 |
| SOFR-OIS-9Y | 26,686.6 | 26,684.5 | 26,686.6 |
| SOFR-OIS-7Y | -17,807.6 | -17,793.6 | -17,807.5 |

The long-dated risk lands on `FF-OIS-30Y` and `FF-SOFR-30Y`, not on a SOFR swap: there is no 30-year SOFR
swap in the curve, so those two instruments together *are* the 30-year SOFR point. A ladder that showed
"30-year SOFR" would be showing an instrument nobody can trade.

## 3. Volatility: the swaption cube

**Product view.** A European swaption is quoted as a normal (Bachelier) vol of the forward swap rate. The
cube holds one SABR smile per (expiry, tenor): α sets the level, ρ the skew, ν the curvature, with β = 0.5
fixed as policy. α is re-solved from the at-the-money quote at every step of the fit, so every node reprices
its ATM vol exactly and only ρ and ν are searched.

- All 81 nodes fit with no parameter on a bound. Worst residual 0.016 bp
  (3Mx30Y); the quotes are rounded to 0.01 bp, so this is the rounding of the input.
- 3–5 solver iterations per node, from a starting point computed from the smile's
  skew and curvature at the money.
- Normal ↔ shifted-lognormal conversion round-trips to 1e-12. Normal → lognormal can fail (a lognormal price
  is bounded, a normal one is not) and that is reported rather than clamped.

**Where the model breaks.** Hagan's formula is an expansion, not a model: far enough below the forward its
implied density goes negative, which is a butterfly arbitrage. 50 of the 81 nodes have a
negative density somewhere inside four standard deviations. The first few, by expiry:

| Node | Forward % | ATM vol (bp) | ρ | ν | Fit residual (bp) | Density < 0 below strike % |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 1Yx1Y | 3.3252 | 98.35 | -0.0991 | 0.4628 | 0.0041 | 0.046 |
| 1Yx2Y | 3.3843 | 95.81 | -0.1034 | 0.4510 | 0.0030 | 0.028 |
| 1Yx3Y | 3.4297 | 93.86 | -0.1072 | 0.4412 | 0.0017 | 0.027 |
| 1Yx5Y | 3.4961 | 91.16 | -0.1139 | 0.4269 | 0.0036 | 0.027 |
| 1Yx7Y | 3.5442 | 89.47 | -0.1193 | 0.4173 | 0.0046 | 0.027 |
| 2Yx1Y | 3.3726 | 92.94 | -0.1332 | 0.3398 | 0.0035 | 0.053 |
| 2Yx2Y | 3.4292 | 90.40 | -0.1380 | 0.3282 | 0.0034 | 0.052 |
| 2Yx3Y | 3.4724 | 88.45 | -0.1415 | 0.3185 | 0.0029 | 0.052 |
| 2Yx5Y | 3.5348 | 85.75 | -0.1485 | 0.3036 | 0.0037 | 0.030 |
| 2Yx7Y | 3.5796 | 84.06 | -0.1538 | 0.2943 | 0.0032 | 0.030 |
| 2Yx10Y | 3.6295 | 82.59 | -0.1606 | 0.2850 | 0.0052 | 0.030 |
| 2Yx15Y | 3.6866 | 81.45 | -0.1675 | 0.2784 | 0.0047 | 0.030 |

![sabr_smile.png](img/sabr_smile.png)

The fit residual says nothing about this: the stressed 10Y×30Y node fits its seven quotes as well as any
other. The check has to be made on the density. It matters for anything that integrates the smile over
strikes (CMS, digitals), and it bounds how far from the money the cube should be read.

## 4. Model: LGM one-factor (Hull–White), calibrated to the co-terminal strip

**Model view.** One Gaussian factor. In the LGM parameterisation the model is two functions: H(t), fixed by
the mean reversion κ, and ζ(t), the variance of the state. A swaption expiring at t depends on ζ only
through ζ(t), so the nine co-terminal swaptions of the 10-year callable swap pin ζ one point at a time by a
1-D root search. The strike is the trade's own (3.50%), read off the smile.

Calibration at κ = 3%:

| Expiry date | Expiry (y) | Tenor (y) | Forward % | Market vol (bp) | Model vol (bp) | Error (bp) | σ (bp) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 2027-09-09 | 1.01 | 9.0 | 3.4216 | 88.67 | 88.67 | 3e-13 | 100.5 |
| 2028-09-11 | 2.02 | 8.0 | 3.4364 | 83.71 | 83.71 | 3e-13 | 89.7 |
| 2029-09-10 | 3.02 | 7.0 | 3.5361 | 80.30 | 80.30 | -5e-13 | 83.1 |
| 2030-09-09 | 4.02 | 6.0 | 3.5842 | 78.61 | 78.61 | -4e-13 | 83.7 |
| 2031-09-09 | 5.02 | 5.0 | 3.6280 | 77.86 | 77.86 | 2e-13 | 85.1 |
| 2032-09-09 | 6.02 | 4.0 | 3.6643 | 78.16 | 78.16 | -9e-14 | 89.8 |
| 2033-09-09 | 7.02 | 3.0 | 3.7005 | 78.73 | 78.73 | 6e-13 | 91.9 |
| 2034-09-11 | 8.02 | 2.0 | 3.7378 | 80.22 | 80.22 | 4e-13 | 99.1 |
| 2035-09-10 | 9.02 | 1.0 | 3.7739 | 82.38 | 82.38 | 3e-13 | 106.3 |

Worst repricing error over κ = 0%, 3%, 10%: 1e-12 bp.

![lgm_sigma.png](img/lgm_sigma.png)

*The strip does not identify κ.* All three models reprice the same nine swaptions; they differ in how much
short-rate volatility is needed to do it. Higher κ makes long-dated rates less sensitive to the state, so σ
must rise to keep the same swaption vols. What κ changes is the correlation of rates *across* exercise
dates, which the Europeans cannot see and the Bermudan can (section 6).

## 5. Product: the Bermudan swaption, priced three ways

**Product view.** A 10-year fixed-for-SOFR swap that one side may enter on any of nine annual dates from
year 1 ("10 non-call 1"). On each date the holder compares the value of entering the remaining swap with
the value of keeping the option. Prices are in basis points of notional; κ = 3%.

| Trade | Best European | Lower bound (LSMC) | PDE | Upper bound (A–B) | Duality gap | P(exercise) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 10nc1 payer 2.50% (ITM) | 738.8 | 851.60 ± 0.14 | 851.42 | 852.33 ± 0.17 | 0.73 ± 0.08 | 91% |
| 10nc1 payer 3.50% (near ATM) | 322.8 | 442.91 ± 0.11 | 442.88 | 443.35 ± 0.12 | 0.43 ± 0.05 | 76% |
| 10nc1 payer 4.50% (OTM) | 144.6 | 222.99 ± 0.10 | 222.98 | 223.17 ± 0.11 | 0.18 ± 0.03 | 59% |
| 10nc1 receiver 3.50% | 331.4 | 431.74 ± 0.16 | 431.91 | 432.29 ± 0.18 | 0.55 ± 0.07 | 61% |

- **Lower bound.** Longstaff–Schwartz: fit an exercise rule by regression on 20,000 paths, then value that
  rule on 4 million independent paths. Any rule is a strategy the holder could follow, so it cannot be
  worth more than the option.
- **Upper bound.** Andersen–Broadie: the dual formulation, with the martingale built from the same rule by
  nested simulation (2,000 outer × 500 inner paths).
- **PDE.** Crank–Nicolson on a 1600 × 1600 grid, converged to about 0.003 bp, and confirmed to 0.001 bp
  by an independent backward induction that integrates against the exact Gaussian transition density.
- The three agree on every trade: the PDE lies within two standard errors of both bounds (the largest
  excursion is 1.3 standard errors), and the duality gap, which is measured much more precisely
  than either bound, is
  0.2–0.7 bp, about 0.1% of the price.
- The Bermudan is worth far more than its best single European: the right to choose the date is a large
  part of the value, more so out of the money.

Wall-clock per pricing (PDE at 400 × 400, which is already within 0.05 bp):

| Trade | PDE (ms) | LSMC, 4M paths (ms) | Upper bound (ms) |
| --- | ---: | ---: | ---: |
| 10nc1 payer 2.50% (ITM) | 2.3 | 1745 | 1613 |
| 10nc1 payer 3.50% (near ATM) | 2.3 | 2446 | 2119 |
| 10nc1 payer 4.50% (OTM) | 2.3 | 2645 | 2871 |
| 10nc1 receiver 3.50% | 2.3 | 2491 | 2437 |

### 5.1 PDE: what converges, and what Rannacher start-up actually buys

Under LGM the deflated value satisfies the heat equation in the variance clock ζ, with no drift or
discounting terms: the change of numeraire has absorbed them.

| Steps (space = time) | Error with Rannacher (bp) | Error, Crank–Nicolson only (bp) | Time (ms) |
| --- | ---: | ---: | ---: |
| 50 | 4.9e+00 | 4.4e+00 | 0.07 |
| 100 | 4.4e-01 | 3.1e-01 | 0.21 |
| 200 | 4.1e-01 | 3.8e-01 | 0.67 |
| 400 | 4.8e-02 | 3.8e-02 | 2.28 |
| 800 | 5.0e-03 | 2.7e-03 | 8.47 |
| 1600 | 3.7e-03 | 3.1e-03 | 33.10 |
| 3200 | 5.0e-04 | 3.5e-04 | 134.06 |

![pde_convergence.png](img/pde_convergence.png)

- **Second order overall, not monotone.** The error falls 10,000-fold over a 64-fold refinement, but
  unevenly: it depends on where the exercise boundary falls between grid nodes, so one doubling can do
  nothing and the next cut the error sixteen-fold. Convergence has to be judged over several doublings.
- **Rannacher start-up costs price accuracy here.** The textbook advice is to replace the first steps after
  each exercise date with implicit ones, to damp the oscillations Crank–Nicolson would otherwise carry. For
  this payoff the price is *more* accurate without it (middle panel): exercising creates a kink, not a
  jump, and the implicit steps bring their own time-stepping error.
- **It buys smooth greeks.** With few time steps on a fine space grid, the second derivative from plain
  Crank–Nicolson oscillates; with the start-up it does not (right panel). It is on by default because the
  output is hedged, not just quoted.

### 5.2 Monte Carlo: the exercise rule, the bias, the error

PDE reference 442.88 bp.

| Regressors | Degree | Lower bound (bp) | Duality gap (bp) |
| --- | ---: | ---: | ---: |
| state | 0 | 433.55 ± 0.36 | 9.90 ± 0.44 |
| state | 1 | 441.94 ± 0.36 | 0.76 ± 0.08 |
| state | 2 | 442.22 ± 0.36 | 0.49 ± 0.06 |
| state | 3 | 442.30 ± 0.35 | 0.39 ± 0.05 |
| state | 4 | 442.31 ± 0.35 | 0.42 ± 0.06 |
| state | 6 | 442.31 ± 0.35 | 0.47 ± 0.06 |
| exercise value | 0 | 433.55 ± 0.36 | 9.90 ± 0.44 |
| exercise value | 1 | 441.23 ± 0.36 | 1.46 ± 0.13 |
| exercise value | 2 | 442.27 ± 0.36 | 0.41 ± 0.06 |
| exercise value | 3 | 442.31 ± 0.35 | 0.43 ± 0.05 |
| exercise value | 4 | 442.30 ± 0.35 | 0.41 ± 0.06 |
| exercise value | 6 | 442.33 ± 0.35 | 0.44 ± 0.05 |

![lsmc.png](img/lsmc.png)

- **Degree 0** (continuation value = a constant, so the rule cannot see the state) loses about 9 bp, and
  the gap shows it without needing the PDE: the gap is a model-free measure of how good the rule is.
- **Degree 1–2 is already enough** in a one-factor model; beyond degree 3 nothing improves. Regressing on
  the swap's value works as well as regressing on the state, and is the choice that carries over to models
  with more factors.
- **In-sample bias.** Valuing the rule on its own training paths overstates the price, by about
  15 bp with 250 paths; on fresh paths
  the same rule is *below* the true price, as a lower bound must be. The two meet as the training set grows.
  (Averages over 200 seeds.)
- **Error** falls as 1/√N; antithetic paths cut it by 1.8×, the equivalent of 3.1× the paths.

## 6. What drives the Bermudan: model and product behaviour

![bermudan_behaviour.png](img/bermudan_behaviour.png)

**Mean reversion.** The calibration swaptions do not move with κ; the Bermudan does.

| κ | Bermudan (bp) | Best European (bp) | Switch value (bp) |
| --- | ---: | ---: | ---: |
| -2.0% | 420.5 | 322.8 | 97.6 |
| -1.0% | 424.9 | 322.8 | 102.1 |
| 0.0% | 429.4 | 322.8 | 106.5 |
| 1.0% | 433.9 | 322.8 | 111.0 |
| 2.0% | 438.4 | 322.8 | 115.5 |
| 3.0% | 442.9 | 322.8 | 120.0 |
| 5.0% | 451.9 | 322.8 | 129.0 |
| 7.5% | 463.1 | 322.8 | 140.3 |
| 10.0% | 474.4 | 322.8 | 151.6 |
| 15.0% | 496.7 | 322.8 | 173.8 |

From κ = -2% to 15% the switch value goes from 98 to
174 bp with every calibration instrument repriced exactly throughout. Higher κ
means a shock today has decayed by the later exercise dates, so the exercise decisions are less correlated
and the right to choose among them is worth more. κ is the number a Bermudan desk marks, and this table is
its price.

**Exercise rights.** Value is non-decreasing in the set of call dates, by construction of the product; the
middle panel adds the nine dates one at a time from the last.

**Exercise boundary** (payer, strike 3.50%):

| Call date | Forward swap rate % | European (bp) | Exercise if swap rate above % | P(exercise here) |
| --- | ---: | ---: | ---: | ---: |
| 2027-09-09 | 3.4216 | 238.2 | 4.5438 | 11.7% |
| 2028-09-11 | 3.4364 | 289.8 | 4.5935 | 11.2% |
| 2029-09-10 | 3.5361 | 322.8 | 4.4224 | 12.3% |
| 2030-09-09 | 3.5842 | 317.6 | 4.3413 | 8.4% |
| 2031-09-09 | 3.6280 | 294.7 | 4.2527 | 7.7% |
| 2032-09-09 | 3.6643 | 258.2 | 4.1414 | 6.9% |
| 2033-09-09 | 3.7005 | 209.4 | 4.0076 | 6.1% |
| 2034-09-11 | 3.7378 | 150.4 | 3.8103 | 6.2% |
| 2035-09-10 | 3.7739 | 81.1 | 3.5000 | 5.8% |

At the last date the boundary is the strike: nothing is left to wait for. Earlier, the swap must be in the
money by a margin, about 100 bp at the first call date, because exercising gives up the remaining
optionality.

## 7. Bermudan risk

Notional 1,000,000, payer 3.50%, κ = 3%. Each number is a full revaluation: bump the market, re-bootstrap
the curve, **re-calibrate the model**, reprice by PDE. Re-calibrating is the point: the hedge instruments
are the co-terminal swaptions, and the model has to keep repricing them after the bump.

![bermudan_risk.png](img/bermudan_risk.png)

**Delta.** Total 285 per bp against 748 for the swap from the first call date: the option's delta
is about 38% of the underlying's.

| Quote | Bermudan | Forward swap |
| --- | ---: | ---: |
| SR3-U6 | -4.9 | -25.0 |
| SR3-Z6 | -4.9 | -25.0 |
| SR3-H7 | -4.9 | -25.0 |
| SR3-M7 | -4.7 | -23.3 |
| SOFR-OIS-3Y | -31.2 | 0.2 |
| SOFR-OIS-4Y | -31.7 | 0.2 |
| SOFR-OIS-5Y | -30.8 | 0.1 |
| SOFR-OIS-6Y | -32.8 | 0.3 |
| SOFR-OIS-7Y | -32.6 | 0.5 |
| SOFR-OIS-8Y | -37.2 | -0.1 |
| SOFR-OIS-9Y | -42.5 | 2.7 |
| SOFR-OIS-10Y | 552.9 | 844.2 |

The swap's delta is one long bucket at 10 years and a short at its 1-year start. The Bermudan has the same
long end and a *short in every intermediate year*: it is a probability-weighted mix of swaps starting in
years 1 to 9, each of which is long the 10-year rate and short the rate to its own start date. No risk
appears beyond 10 years or on the Fed Funds basis.

**Vega**, by calibration instrument:

| Expiry | Bermudan vega | European's own vega | Ratio |
| --- | ---: | ---: | ---: |
| 2027-09-09 | 51.2 | 299.5 | 17% |
| 2028-09-11 | 59.9 | 370.0 | 16% |
| 2029-09-10 | 89.0 | 389.1 | 23% |
| 2030-09-09 | 83.9 | 377.6 | 22% |
| 2031-09-09 | 75.0 | 344.8 | 22% |
| 2032-09-09 | 63.7 | 296.2 | 22% |
| 2033-09-09 | 54.8 | 235.3 | 23% |
| 2034-09-11 | 44.2 | 164.1 | 27% |
| 2035-09-10 | 28.4 | 85.4 | 33% |

Total 550 against 2562 for holding all nine Europeans. Each ratio is roughly the hedge ratio in that
European; they are all well below one because the Bermudan is one option on the best of nine, not nine
options.

## 8. Exposure: the same model, seen from counterparty risk

The expected positive exposure of a swap at a future date is E[max(V, 0)], which is by definition the price
of a European swaption on the remaining swap. A swap's exposure profile is a strip of co-terminal
swaptions: the same instruments the Bermudan is calibrated to. So with this calibration the profile below
is not a model output at all from year 1 on. It is the market's swaption prices.

| Years | EE, closed form (bp) | EE, Monte Carlo (bp) | PFE 97.5% (bp) | ENE (bp) |
| --- | ---: | ---: | ---: | ---: |
| 0.01 | 18.0 | 18.0 ± 0.2 | 130 | -60.3 |
| 1.01 | 238.2 | 237.2 ± 1.6 | 1168 | -296.9 |
| 2.02 | 289.8 | 286.2 ± 1.8 | 1334 | -331.4 |
| 3.02 | 322.8 | 325.4 ± 1.9 | 1351 | -302.5 |
| 4.02 | 317.6 | 316.4 ± 1.7 | 1255 | -277.8 |
| 5.02 | 294.7 | 291.1 ± 1.6 | 1103 | -245.2 |
| 6.02 | 258.2 | 259.1 ± 1.3 | 915 | -208.3 |
| 7.02 | 209.4 | 209.5 ± 1.0 | 701 | -164.5 |
| 8.03 | 150.4 | 150.7 ± 0.7 | 473 | -115.7 |
| 9.02 | 81.1 | 81.6 ± 0.4 | 237 | -61.4 |

![exposure.png](img/exposure.png)

- The closed form and a 50,000-path simulation agree within sampling error at every date, which is how an
  exposure engine can be validated against a pricing model.
- The hump is diffusion against amortisation: rate uncertainty grows with time while the remaining annuity
  shrinks.
- CVA on this swap is this strip weighted by default probabilities. The right to cancel the swap is the
  Bermudan of section 5 on the offsetting swap. The exposure world and the exotics world are pricing the
  same options.

Not modelled: collateral, netting, margin period of risk, wrong-way risk.
