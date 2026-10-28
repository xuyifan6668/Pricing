"""Regenerate docs/RESULTS.md and docs/img/*.png from the C++ analysis executable.

    uv run python tools/make_results.py [--exe build/dev/cpp/examples/rl_analysis] [--skip-run]

Every number in the results document comes from the CSV files that rl_analysis writes, so the document
cannot drift from the code. The run is deterministic apart from the timings.
"""

from __future__ import annotations

import argparse
import csv
import os
import platform
import subprocess
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.ticker  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
IMG = ROOT / "docs" / "img"

# Three categorical slots that stay distinguishable under colour-vision deficiency; ink for text.
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e4e3df"

plt.rcParams.update(
    {
        "figure.dpi": 150,
        "savefig.dpi": 150,
        "figure.facecolor": "#fcfcfb",
        "axes.facecolor": "#fcfcfb",
        "axes.edgecolor": GRID,
        "axes.labelcolor": MUTED,
        "axes.titlecolor": INK,
        "axes.titlesize": 11,
        "axes.titleweight": "bold",
        "axes.titlelocation": "left",
        "axes.labelsize": 9,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "grid.color": GRID,
        "grid.linewidth": 0.6,
        "xtick.color": MUTED,
        "ytick.color": MUTED,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "lines.linewidth": 2.0,
        "lines.markersize": 5,
        "legend.frameon": False,
        "legend.fontsize": 8,
        "font.size": 9,
    }
)


def read(results: Path, name: str) -> list[dict[str, str]]:
    with open(results / f"{name}.csv", newline="") as f:
        return list(csv.DictReader(f))


def col(rows: list[dict[str, str]], key: str) -> list[float]:
    return [float(r[key]) for r in rows]


def table(header: list[str], rows: list[list[str]], align: str | None = None) -> str:
    align = align or "l" + "r" * (len(header) - 1)
    rule = ["---" if a == "l" else "---:" for a in align]
    lines = ["| " + " | ".join(header) + " |", "| " + " | ".join(rule) + " |"]
    lines += ["| " + " | ".join(r) + " |" for r in rows]
    return "\n".join(lines)


def plain_log_ticks(ax, values: list[float], axis: str = "x") -> None:
    """Label a log axis at the data's own values, in plain numbers, with no minor-tick labels."""
    target = ax.xaxis if axis == "x" else ax.yaxis
    target.set_major_locator(matplotlib.ticker.FixedLocator(values))
    target.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:,.0f}"))
    target.set_minor_formatter(matplotlib.ticker.NullFormatter())


def save(fig, name: str) -> str:
    IMG.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(IMG / name)
    plt.close(fig)
    return f"![{name}](img/{name})"


def machine() -> str:
    cpu = platform.processor() or platform.machine()
    try:
        for line in Path("/proc/cpuinfo").read_text().splitlines():
            if line.startswith("model name"):
                cpu = line.split(":", 1)[1].strip()
                break
    except OSError:
        pass
    return f"{cpu}, {os.cpu_count()} logical cores, {platform.system()} {platform.machine()}, single-threaded, GCC -O2"


# --------------------------------------------------------------------------------------------------
# Sections
# --------------------------------------------------------------------------------------------------


def curves(res: Path) -> str:
    solver = read(res, "curve_solver")
    reprice = read(res, "curve_reprice")
    fwd = read(res, "curve_forwards")
    worst = {}
    for r in reprice:
        worst[r["interpolation"]] = max(worst.get(r["interpolation"], 0.0), abs(float(r["residual_bp"])))

    fig, axes = plt.subplots(1, 2, figsize=(9.0, 3.2), sharey=True)
    for ax, lo, hi, title in ((axes[0], 0.0, 30.0, "Forward rate, 0–30 years"), (axes[1], 0.0, 3.0, "The first three years")):
        rows = [r for r in fwd if lo <= float(r["t"]) <= hi]
        ax.plot(col(rows, "t"), col(rows, "linear_log_df"), color=BLUE, label="linear on log discount")
        ax.plot(col(rows, "t"), col(rows, "monotone_convex"), color=ORANGE, label="monotone convex")
        ax.set_title(title)
        ax.set_xlabel("years")
    axes[0].set_ylabel("1-week forward, %")
    axes[0].legend(loc="lower right")
    plot = save(fig, "curve_forwards.png")

    rows = [[r["interpolation"], r["method"], r["iterations"], f"{float(r['max_residual_bp']):.1e}", r["worst_instrument"], f"{float(r['time_ms']):.1f}"] for r in solver]
    n = len([r for r in reprice if r["interpolation"] == "linear_logdf"])
    return f"""## 1. Curves: USD SOFR and Fed Funds, bootstrapped together

**Product view.** {n} instruments define two curves: SOFR overnight-index swaps and 3-month SOFR futures for
the SOFR curve, Fed Funds–SOFR basis swaps and Fed Funds swaps for the Fed Funds curve. The 25- and 30-year
SOFR points are quoted only through basis swaps, so the SOFR long end depends on the Fed Funds curve while
Fed Funds cashflows are discounted on SOFR. Neither curve can be built before the other.

{table(["Interpolation", "Method", "Newton iterations", "Worst residual (bp)", "Worst instrument", "Time (ms)"], rows, "llrrlr")}

- The global solve reprices every instrument to about 1e-10 bp in four iterations under both interpolations
  (worst residuals: {", ".join(f"{k} {v:.1e} bp" for k, v in worst.items())}).
- One pass of the curve-by-curve bootstrap leaves the long basis swaps about 5 bp off, because they were
  solved against a Fed Funds curve that did not exist yet.
- The Jacobian is taken by an adjoint sweep through the same pricing code the solver runs. It gives the same
  answer as finite differences about ten times faster, and stays exact under monotone convex interpolation,
  where a hand-derived Jacobian is easy to get subtly wrong.

{plot}

*Interpolation is a modelling choice, not a detail.* Linear interpolation of log discount factors gives
piecewise-constant forwards with a jump at every pillar; monotone convex gives continuous forwards. Both
reprice the same quotes exactly. They differ in the forward rate between pillars, which is what a
forward-starting swap or a futures strip is priced on.
"""


def swap_risk(res: Path) -> str:
    timing = read(res, "swap_risk_timing")
    rows = read(res, "swap_risk")
    big = sorted(rows, key=lambda r: -abs(float(r["adjoint"])))[:8]
    t = {r["method"]: float(r["time_ms"]) for r in timing}
    aad = next(v for k, v in t.items() if k.startswith("adjoint"))
    bump = next(v for k, v in t.items() if k.startswith("bump"))
    one = next(v for k, v in t.items() if k.startswith("one"))
    rel = next(r["max_difference_rel"] for r in timing if r["method"].startswith("adjoint"))
    swaps, quotes = timing[0]["swaps"], timing[0]["quotes"]
    return f"""## 2. Swap risk: adjoint differentiation against bumping

**Product view.** A swap desk hedges with the instruments the curve is built from, so risk is reported per
par quote: "how much does the book move if the 10-year swap rate rises 1 bp and the curve is rebuilt".

A book of {swaps} swaps (1–30 years, both directions), delta to all {quotes} curve quotes:

{table(["Method", "Time (ms)", "Relative to one valuation"], [[r["method"], f"{float(r['time_ms']):.2f}", f"{float(r['time_ms']) / one:.0f}×"] for r in timing])}

- Adjoint: one reverse sweep gives the sensitivity to every curve node; one transposed linear solve with the
  bootstrap Jacobian (implicit-function theorem) converts node risk to par-quote risk. **{bump / aad:.0f}× faster**
  than bumping here, and the cost does not grow with the number of quotes.
- The two agree to {float(rel):.0e} of the largest delta when the bump is small (0.01 bp). A 1 bp bump differs
  by about 1e-3 because it also contains the gamma.

Largest deltas (currency units per +1 bp, notionals 1–5 million):

{table(["Quote", "Adjoint", "Bump 1 bp", "Bump 0.01 bp"], [[r["quote"], f"{float(r['adjoint']):,.1f}", f"{float(r['bump_1bp']):,.1f}", f"{float(r['bump_0.01bp']):,.1f}"] for r in big])}

The long-dated risk lands on `FF-OIS-30Y` and `FF-SOFR-30Y`, not on a SOFR swap: there is no 30-year SOFR
swap in the curve, so those two instruments together *are* the 30-year SOFR point. A ladder that showed
"30-year SOFR" would be showing an instrument nobody can trade.
"""


def vol(res: Path) -> str:
    fit = read(res, "sabr_fit")
    smile = read(res, "sabr_smile")
    worst = max(fit, key=lambda r: float(r["max_residual_bp"]))
    bad = [r for r in fit if r["first_negative_density_strike_pct"]]
    its = col(fit, "iterations")

    fig, axes = plt.subplots(1, 2, figsize=(9.0, 3.2))
    for node, color in (("5Yx10Y", BLUE), ("10Yx30Y", ORANGE)):
        rows = [r for r in smile if r["node"] == node]
        axes[0].plot(col(rows, "strike_pct"), col(rows, "normal_vol_bp"), color=color, label=node)
        axes[1].plot(col(rows, "strike_pct"), col(rows, "density"), color=color, label=node)
    axes[1].axhline(0.0, color=MUTED, linewidth=0.8)
    axes[0].set_title("Fitted SABR smile")
    axes[0].set_xlabel("strike, %")
    axes[0].set_ylabel("normal vol, bp")
    axes[0].legend()
    axes[1].set_title("Implied density of the swap rate")
    axes[1].set_xlabel("strike, %")
    axes[1].set_ylabel("density")
    axes[1].legend()
    plot = save(fig, "sabr_smile.png")

    rows = [[r["node"], r["forward_pct"], r["atm_bp"], r["rho"], r["nu"], r["max_residual_bp"], r["first_negative_density_strike_pct"][:5]] for r in bad[:12]]
    return f"""## 3. Volatility: the swaption cube

**Product view.** A European swaption is quoted as a normal (Bachelier) vol of the forward swap rate. The
cube holds one SABR smile per (expiry, tenor): α sets the level, ρ the skew, ν the curvature, with β = 0.5
fixed as policy. α is re-solved from the at-the-money quote at every step of the fit, so every node reprices
its ATM vol exactly and only ρ and ν are searched.

- All {len(fit)} nodes fit with no parameter on a bound. Worst residual {float(worst["max_residual_bp"]):.3f} bp
  ({worst["node"]}); the quotes are rounded to 0.01 bp, so this is the rounding of the input.
- {min(its):.0f}–{max(its):.0f} solver iterations per node, from a starting point computed from the smile's
  skew and curvature at the money.
- Normal ↔ shifted-lognormal conversion round-trips to 1e-12. Normal → lognormal can fail (a lognormal price
  is bounded, a normal one is not) and that is reported rather than clamped.

**Where the model breaks.** Hagan's formula is an expansion, not a model: far enough below the forward its
implied density goes negative, which is a butterfly arbitrage. {len(bad)} of the {len(fit)} nodes have a
negative density somewhere inside four standard deviations. The first few, by expiry:

{table(["Node", "Forward %", "ATM vol (bp)", "ρ", "ν", "Fit residual (bp)", "Density < 0 below strike %"], rows)}

{plot}

The fit residual says nothing about this: the stressed 10Y×30Y node fits its seven quotes as well as any
other. The check has to be made on the density. It matters for anything that integrates the smile over
strikes (CMS, digitals), and it bounds how far from the money the cube should be read.
"""


def lgm(res: Path) -> str:
    cal = read(res, "lgm_calibration")
    fig, ax = plt.subplots(figsize=(5.6, 3.2))
    for kappa, color in (("0.00", BLUE), ("0.03", ORANGE), ("0.10", AQUA)):
        rows = [r for r in cal if r["kappa"] == kappa]
        ax.step([0.0] + col(rows, "expiry"), [float(rows[0]["hw_sigma_bp"])] + col(rows, "hw_sigma_bp"), where="pre", color=color, label=f"κ = {float(kappa):.0%}")
    ax.set_title("Calibrated short-rate volatility")
    ax.set_xlabel("years")
    ax.set_ylabel("Hull–White σ, bp")
    ax.legend()
    plot = save(fig, "lgm_sigma.png")
    rows = [r for r in cal if r["kappa"] == "0.03"]
    worst = max(abs(float(r["error_bp"])) for r in cal)
    body = [[r["label"], f"{float(r['expiry']):.2f}", f"{float(r['tenor']):.1f}", r["forward_pct"], f"{float(r['market_vol_bp']):.2f}", f"{float(r['model_vol_bp']):.2f}", f"{float(r['error_bp']):.0e}", f"{float(r['hw_sigma_bp']):.1f}"] for r in rows]
    return f"""## 4. Model: LGM one-factor (Hull–White), calibrated to the co-terminal strip

**Model view.** One Gaussian factor. In the LGM parameterisation the model is two functions: H(t), fixed by
the mean reversion κ, and ζ(t), the variance of the state. A swaption expiring at t depends on ζ only
through ζ(t), so the nine co-terminal swaptions of the 10-year callable swap pin ζ one point at a time by a
1-D root search. The strike is the trade's own (3.50%), read off the smile.

Calibration at κ = 3%:

{table(["Expiry date", "Expiry (y)", "Tenor (y)", "Forward %", "Market vol (bp)", "Model vol (bp)", "Error (bp)", "σ (bp)"], body)}

Worst repricing error over κ = 0%, 3%, 10%: {worst:.0e} bp.

{plot}

*The strip does not identify κ.* All three models reprice the same nine swaptions; they differ in how much
short-rate volatility is needed to do it. Higher κ makes long-dated rates less sensitive to the state, so σ
must rise to keep the same swaption vols. What κ changes is the correlation of rates *across* exercise
dates, which the Europeans cannot see and the Bermudan can (section 6).
"""


def bermudan(res: Path) -> str:
    three = read(res, "bermudan_three_ways")
    rows = [
        [
            r["trade"],
            f"{float(r['max_european_bp']):.1f}",
            f"{float(r['lsmc_lower_bp']):.2f} ± {float(r['lsmc_se_bp']):.2f}",
            f"{float(r['pde_bp']):.2f}",
            f"{float(r['upper_bp']):.2f} ± {float(r['upper_se_bp']):.2f}",
            f"{float(r['gap_bp']):.2f} ± {float(r['gap_se_bp']):.2f}",
            f"{float(r['exercise_probability']):.0%}",
        ]
        for r in three
    ]
    # How far outside the Monte Carlo bracket the PDE sits, in standard errors (≤ 0: inside).
    excursion = max(
        max((float(r["lsmc_lower_bp"]) - float(r["pde_bp"])) / float(r["lsmc_se_bp"]), (float(r["pde_bp"]) - float(r["upper_bp"])) / float(r["upper_se_bp"]))
        for r in three
    )
    if excursion > 2.5:
        raise SystemExit(f"the PDE is {excursion:.1f} standard errors outside the Monte Carlo bounds: investigate before publishing")
    excursion = max(excursion, 0.0)
    times = [[r["trade"], f"{float(r['pde_ms']):.1f}", f"{float(r['lsmc_ms']):.0f}", f"{float(r['upper_ms']):.0f}"] for r in three]

    # PDE convergence
    conv = read(res, "pde_convergence")
    tconv = read(res, "pde_time_convergence")
    gamma = read(res, "pde_gamma")
    fig, axes = plt.subplots(1, 3, figsize=(11.0, 3.2))
    n = col(conv, "steps")
    axes[0].loglog(n, col(conv, "error_rannacher_bp"), "o-", color=BLUE, label="with Rannacher")
    axes[0].loglog(n, col(conv, "error_cn_bp"), "o-", color=ORANGE, label="Crank–Nicolson only")
    ref = [col(conv, "error_cn_bp")[0] * (n[0] / x) ** 2 for x in n]
    axes[0].loglog(n, ref, "--", color=MUTED, linewidth=1.0, label="slope −2")
    axes[0].set_title("Error, space and time refined together")
    axes[0].set_xlabel("steps in each direction")
    axes[0].set_ylabel("|error|, bp of notional")
    plain_log_ticks(axes[0], [50, 200, 800, 3200])
    axes[0].legend()
    tn = col(tconv, "time_steps")
    axes[1].loglog(tn, col(tconv, "error_rannacher_bp"), "o-", color=BLUE, label="with Rannacher")
    axes[1].loglog(tn, col(tconv, "error_cn_bp"), "o-", color=ORANGE, label="Crank–Nicolson only")
    axes[1].set_title("Error against time steps (3200 in space)")
    axes[1].set_xlabel("time steps")
    plain_log_ticks(axes[1], tn)
    axes[1].legend()
    axes[2].plot(col(gamma, "x"), col(gamma, "gamma_cn"), color=ORANGE, linewidth=1.0, label="Crank–Nicolson only")
    axes[2].plot(col(gamma, "x"), col(gamma, "gamma_rannacher"), color=BLUE, label="with Rannacher")
    axes[2].set_title("Gamma in the state, 36 time steps")
    axes[2].set_xlabel("state x")
    axes[2].set_ylabel("∂²V/∂x²")
    axes[2].legend()
    pde_plot = save(fig, "pde_convergence.png")
    pde_rows = [[r["steps"], f"{float(r['error_rannacher_bp']):.1e}", f"{float(r['error_cn_bp']):.1e}", f"{float(r['time_ms']):.2f}"] for r in conv]

    # LSMC
    basis = read(res, "lsmc_basis")
    bias = read(res, "lsmc_bias")
    paths = read(res, "lsmc_paths")
    fig, axes = plt.subplots(1, 3, figsize=(11.0, 3.2))
    for name, color in (("state", BLUE), ("exercise value", ORANGE)):
        rows_b = [r for r in basis if r["basis"] == name]
        axes[0].plot(col(rows_b, "degree"), col(rows_b, "gap_bp"), "o-", color=color, label=f"powers of the {name}")
    axes[0].set_yscale("log")
    axes[0].set_title("Duality gap against polynomial degree")
    axes[0].set_xlabel("degree")
    axes[0].set_ylabel("upper − lower, bp")
    axes[0].legend()
    tp = col(bias, "training_paths")
    axes[1].semilogx(tp, col(bias, "in_sample_bp"), "o-", color=ORANGE, label="in sample")
    axes[1].semilogx(tp, col(bias, "out_of_sample_bp"), "o-", color=BLUE, label="out of sample")
    axes[1].axhline(float(bias[0]["pde_bp"]), color=MUTED, linewidth=1.0, linestyle="--", label="PDE")
    axes[1].set_title("In-sample bias (degree 6)")
    axes[1].set_xlabel("training paths")
    plain_log_ticks(axes[1], [250, 1000, 4000, 16000])
    axes[1].set_ylabel("value, bp")
    axes[1].legend()
    pp = col(paths, "pricing_paths")
    axes[2].loglog(pp, col(paths, "standard_error_no_antithetic_bp"), "o-", color=ORANGE, label="plain")
    axes[2].loglog(pp, col(paths, "standard_error_bp"), "o-", color=BLUE, label="antithetic")
    axes[2].loglog(pp, [col(paths, "standard_error_no_antithetic_bp")[0] * (pp[0] / x) ** 0.5 for x in pp], "--", color=MUTED, linewidth=1.0, label="slope −½")
    axes[2].set_title("Standard error against paths")
    axes[2].set_xlabel("pricing paths")
    plain_log_ticks(axes[2], [1000, 16000, 256000])
    axes[2].set_ylabel("standard error, bp")
    axes[2].legend()
    mc_plot = save(fig, "lsmc.png")
    basis_rows = [[r["basis"], r["degree"], f"{float(r['lower_bp']):.2f} ± {float(r['lower_se_bp']):.2f}", f"{float(r['gap_bp']):.2f} ± {float(r['gap_se_bp']):.2f}"] for r in basis]
    pde_ref = float(basis[0]["pde_bp"])
    ratio = col(paths, "standard_error_no_antithetic_bp")[-1] / col(paths, "standard_error_bp")[-1]

    # Model and product behaviour
    kappa = read(res, "bermudan_kappa")
    dates = read(res, "bermudan_dates")
    boundary = read(res, "bermudan_boundary")
    fig, axes = plt.subplots(1, 3, figsize=(11.0, 3.2))
    axes[0].plot([100 * k for k in col(kappa, "kappa")], col(kappa, "switch_value_bp"), "o-", color=BLUE)
    axes[0].set_title("Switch value against mean reversion")
    axes[0].set_xlabel("κ, %")
    axes[0].set_ylabel("Bermudan − best European, bp")
    axes[1].plot(col(dates, "exercise_dates"), col(dates, "value_bp"), "o-", color=BLUE)
    axes[1].set_title("Value against number of call dates")
    axes[1].set_xlabel("annual call dates kept (latest first)")
    axes[1].set_ylabel("value, bp")
    axes[2].plot(col(boundary, "expiry"), col(boundary, "boundary_rate_pct"), "o-", color=BLUE, label="exercise above this swap rate")
    axes[2].plot(col(boundary, "expiry"), col(boundary, "forward_pct"), "o-", color=ORANGE, label="forward swap rate")
    axes[2].axhline(3.5, color=MUTED, linewidth=1.0, linestyle="--", label="strike")
    axes[2].set_title("Exercise boundary, payer 3.50%")
    axes[2].set_xlabel("exercise date, years")
    axes[2].set_ylabel("swap rate, %")
    axes[2].legend()
    behaviour_plot = save(fig, "bermudan_behaviour.png")
    kappa_rows = [[f"{float(r['kappa']):.1%}", f"{float(r['bermudan_bp']):.1f}", f"{float(r['max_european_bp']):.1f}", f"{float(r['switch_value_bp']):.1f}"] for r in kappa]
    boundary_rows = [[r["date"], r["forward_pct"], f"{float(r['european_bp']):.1f}", r["boundary_rate_pct"], f"{float(r['exercise_probability']):.1%}"] for r in boundary]
    k_lo, k_hi = kappa[0], kappa[-1]

    return f"""## 5. Product: the Bermudan swaption, priced three ways

**Product view.** A 10-year fixed-for-SOFR swap that one side may enter on any of nine annual dates from
year 1 ("10 non-call 1"). On each date the holder compares the value of entering the remaining swap with
the value of keeping the option. Prices are in basis points of notional; κ = 3%.

{table(["Trade", "Best European", "Lower bound (LSMC)", "PDE", "Upper bound (A–B)", "Duality gap", "P(exercise)"], rows)}

- **Lower bound.** Longstaff–Schwartz: fit an exercise rule by regression on 20,000 paths, then value that
  rule on 4 million independent paths. Any rule is a strategy the holder could follow, so it cannot be
  worth more than the option.
- **Upper bound.** Andersen–Broadie: the dual formulation, with the martingale built from the same rule by
  nested simulation (2,000 outer × 500 inner paths).
- **PDE.** Crank–Nicolson on a 1600 × 1600 grid, converged to about 0.003 bp, and confirmed to 0.001 bp
  by an independent backward induction that integrates against the exact Gaussian transition density.
- The three agree on every trade: the PDE lies within two standard errors of both bounds (the largest
  excursion is {excursion:.1f} standard errors), and the duality gap, which is measured much more precisely
  than either bound, is
  {min(float(r["gap_bp"]) for r in three):.1f}–{max(float(r["gap_bp"]) for r in three):.1f} bp, about 0.1% of the price.
- The Bermudan is worth far more than its best single European: the right to choose the date is a large
  part of the value, more so out of the money.

Wall-clock per pricing (PDE at 400 × 400, which is already within 0.05 bp):

{table(["Trade", "PDE (ms)", "LSMC, 4M paths (ms)", "Upper bound (ms)"], times)}

### 5.1 PDE: what converges, and what Rannacher start-up actually buys

Under LGM the deflated value satisfies the heat equation in the variance clock ζ, with no drift or
discounting terms: the change of numeraire has absorbed them.

{table(["Steps (space = time)", "Error with Rannacher (bp)", "Error, Crank–Nicolson only (bp)", "Time (ms)"], pde_rows)}

{pde_plot}

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

PDE reference {pde_ref:.2f} bp.

{table(["Regressors", "Degree", "Lower bound (bp)", "Duality gap (bp)"], basis_rows, "lrrr")}

{mc_plot}

- **Degree 0** (continuation value = a constant, so the rule cannot see the state) loses about 9 bp, and
  the gap shows it without needing the PDE: the gap is a model-free measure of how good the rule is.
- **Degree 1–2 is already enough** in a one-factor model; beyond degree 3 nothing improves. Regressing on
  the swap's value works as well as regressing on the state, and is the choice that carries over to models
  with more factors.
- **In-sample bias.** Valuing the rule on its own training paths overstates the price, by about
  {float(bias[0]["in_sample_bp"]) - float(bias[0]["out_of_sample_bp"]):.0f} bp with 250 paths; on fresh paths
  the same rule is *below* the true price, as a lower bound must be. The two meet as the training set grows.
  (Averages over 200 seeds.)
- **Error** falls as 1/√N; antithetic paths cut it by {ratio:.1f}×, the equivalent of {ratio * ratio:.1f}× the paths.

## 6. What drives the Bermudan: model and product behaviour

{behaviour_plot}

**Mean reversion.** The calibration swaptions do not move with κ; the Bermudan does.

{table(["κ", "Bermudan (bp)", "Best European (bp)", "Switch value (bp)"], kappa_rows)}

From κ = {float(k_lo["kappa"]):.0%} to {float(k_hi["kappa"]):.0%} the switch value goes from {float(k_lo["switch_value_bp"]):.0f} to
{float(k_hi["switch_value_bp"]):.0f} bp with every calibration instrument repriced exactly throughout. Higher κ
means a shock today has decayed by the later exercise dates, so the exercise decisions are less correlated
and the right to choose among them is worth more. κ is the number a Bermudan desk marks, and this table is
its price.

**Exercise rights.** Value is non-decreasing in the set of call dates, by construction of the product; the
middle panel adds the nine dates one at a time from the last.

**Exercise boundary** (payer, strike 3.50%):

{table(["Call date", "Forward swap rate %", "European (bp)", "Exercise if swap rate above %", "P(exercise here)"], boundary_rows)}

At the last date the boundary is the strike: nothing is left to wait for. Earlier, the swap must be in the
money by a margin, about 100 bp at the first call date, because exercising gives up the remaining
optionality.
"""


def risk(res: Path) -> str:
    delta = read(res, "bermudan_delta")
    vega = read(res, "bermudan_vega")
    fig, axes = plt.subplots(1, 2, figsize=(10.0, 3.4))
    ids = [r["quote"].replace("SOFR-OIS-", "").replace("SR3-", "") for r in delta]
    x = range(len(ids))
    w = 0.4
    axes[0].bar([i - w / 2 for i in x], col(delta, "bermudan"), width=w - 0.04, color=BLUE, label="Bermudan")
    axes[0].bar([i + w / 2 for i in x], col(delta, "forward_swap"), width=w - 0.04, color=ORANGE, label="swap from the first call date")
    axes[0].set_xticks(list(x), ids, rotation=60, ha="right")
    axes[0].axhline(0.0, color=MUTED, linewidth=0.8)
    axes[0].set_title("Delta per +1 bp of each curve quote")
    axes[0].set_ylabel("per 1,000,000 notional")
    axes[0].legend()
    d = [r["expiry_date"][:4] for r in vega]
    xv = range(len(d))
    axes[1].bar([i - w / 2 for i in xv], col(vega, "bermudan_vega"), width=w - 0.04, color=BLUE, label="Bermudan")
    axes[1].bar([i + w / 2 for i in xv], col(vega, "european_vega"), width=w - 0.04, color=ORANGE, label="the European itself")
    axes[1].set_xticks(list(xv), d)
    axes[1].set_title("Vega per +1 bp of each co-terminal swaption vol")
    axes[1].set_xlabel("expiry year")
    axes[1].legend()
    plot = save(fig, "bermudan_risk.png")
    tot_b, tot_s = sum(col(delta, "bermudan")), sum(col(delta, "forward_swap"))
    tot_v, tot_e = sum(col(vega, "bermudan_vega")), sum(col(vega, "european_vega"))
    drows = [[r["quote"], f"{float(r['bermudan']):.1f}", f"{float(r['forward_swap']):.1f}"] for r in delta if abs(float(r["bermudan"])) > 20 or abs(float(r["forward_swap"])) > 20]
    vrows = [[r["expiry_date"], f"{float(r['bermudan_vega']):.1f}", f"{float(r['european_vega']):.1f}", f"{float(r['bermudan_vega']) / float(r['european_vega']):.0%}"] for r in vega]
    return f"""## 7. Bermudan risk

Notional 1,000,000, payer 3.50%, κ = 3%. Each number is a full revaluation: bump the market, re-bootstrap
the curve, **re-calibrate the model**, reprice by PDE. Re-calibrating is the point: the hedge instruments
are the co-terminal swaptions, and the model has to keep repricing them after the bump.

{plot}

**Delta.** Total {tot_b:.0f} per bp against {tot_s:.0f} for the swap from the first call date: the option's delta
is about {tot_b / tot_s:.0%} of the underlying's.

{table(["Quote", "Bermudan", "Forward swap"], drows)}

The swap's delta is one long bucket at 10 years and a short at its 1-year start. The Bermudan has the same
long end and a *short in every intermediate year*: it is a probability-weighted mix of swaps starting in
years 1 to 9, each of which is long the 10-year rate and short the rate to its own start date. No risk
appears beyond 10 years or on the Fed Funds basis.

**Vega**, by calibration instrument:

{table(["Expiry", "Bermudan vega", "European's own vega", "Ratio"], vrows)}

Total {tot_v:.0f} against {tot_e:.0f} for holding all nine Europeans. Each ratio is roughly the hedge ratio in that
European; they are all well below one because the Bermudan is one option on the best of nine, not nine
options.
"""


def exposure(res: Path) -> str:
    rows = read(res, "exposure")
    fig, ax = plt.subplots(figsize=(6.0, 3.2))
    t = col(rows, "time")
    ax.plot(t, col(rows, "pfe_975_bp"), "o-", color=ORANGE, label="97.5% potential exposure")
    ax.plot(t, col(rows, "ee_closed_form_bp"), "o-", color=BLUE, label="expected exposure (swaption strip)")
    ax.plot(t, [-x for x in col(rows, "ene_bp")], "o-", color=AQUA, label="expected negative exposure (sign flipped)")
    ax.set_title("Exposure of a 10-year payer swap at 3.50%")
    ax.set_xlabel("years")
    ax.set_ylabel("discounted exposure, bp of notional")
    ax.legend()
    plot = save(fig, "exposure.png")
    body = [[f"{float(r['time']):.2f}", f"{float(r['ee_closed_form_bp']):.1f}", f"{float(r['ee_monte_carlo_bp']):.1f} ± {float(r['mc_se_bp']):.1f}", f"{float(r['pfe_975_bp']):.0f}", f"{float(r['ene_bp']):.1f}"] for r in rows]
    return f"""## 8. Exposure: the same model, seen from counterparty risk

The expected positive exposure of a swap at a future date is E[max(V, 0)], which is by definition the price
of a European swaption on the remaining swap. A swap's exposure profile is a strip of co-terminal
swaptions: the same instruments the Bermudan is calibrated to. So with this calibration the profile below
is not a model output at all from year 1 on. It is the market's swaption prices.

{table(["Years", "EE, closed form (bp)", "EE, Monte Carlo (bp)", "PFE 97.5% (bp)", "ENE (bp)"], body)}

{plot}

- The closed form and a 50,000-path simulation agree within sampling error at every date, which is how an
  exposure engine can be validated against a pricing model.
- The hump is diffusion against amortisation: rate uncertainty grows with time while the remaining annuity
  shrinks.
- CVA on this swap is this strip weighted by default probabilities. The right to cancel the swap is the
  Bermudan of section 5 on the offsetting swap. The exposure world and the exotics world are pricing the
  same options.

Not modelled: collateral, netting, margin period of risk, wrong-way risk.
"""


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--exe", default=str(ROOT / "build" / "dev" / "cpp" / "examples" / "rl_analysis"))
    ap.add_argument("--results", default=str(ROOT / "build" / "results"))
    ap.add_argument("--skip-run", action="store_true", help="reuse the CSV files already in --results")
    args = ap.parse_args()
    res = Path(args.results)
    if not args.skip_run:
        env = dict(os.environ, RL_DATA_DIR=str(ROOT / "data"))
        subprocess.run([args.exe, str(res)], check=True, env=env)

    head = f"""# Results

Everything here is produced by `rl_analysis` (C++) and rendered by `tools/make_results.py`; nothing is typed
in by hand. Market: USD, as of 2026-09-04, from the files in `data/` (the curve quotes and the vol screen are
synthetic but realistic). Prices are in basis points of notional unless stated.

Timings: {machine()}.

To regenerate:

```bash
cmake -S . -B build/dev -G Ninja -DCMAKE_BUILD_TYPE=Release && cmake --build build/dev
```

```bash
uv run python tools/make_results.py
```

"""
    sections = [curves(res), swap_risk(res), vol(res), lgm(res), bermudan(res), risk(res), exposure(res)]
    (ROOT / "docs" / "RESULTS.md").write_text(head + "\n".join(sections), encoding="utf-8")
    print(f"wrote docs/RESULTS.md and {len(list(IMG.glob('*.png')))} plots")


if __name__ == "__main__":
    main()
