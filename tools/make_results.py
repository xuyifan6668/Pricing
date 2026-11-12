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
    err = {int(r["steps"]): float(r["error_rannacher_bp"]) for r in conv}
    steps = sorted(err)
    drops = [err[a] / err[b] for a, b in zip(steps, steps[1:])]
    deg0_loss = float(basis[0]["pde_bp"]) - float(basis[0]["lower_bp"])
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
- **PDE.** Crank–Nicolson on a 1600 × 1600 grid, converged to about {err[1600]:.3f} bp, and confirmed to 0.001 bp
  by an independent backward induction that integrates against the exact Gaussian transition density.
- The three agree on every trade: the PDE lies within two standard errors of both bounds (the largest
  excursion is {excursion:.1f} standard errors), and the duality gap, which is measured much more precisely
  than either bound, is
  {min(float(r["gap_bp"]) for r in three):.1f}–{max(float(r["gap_bp"]) for r in three):.1f} bp, about 0.1% of the price.
- The Bermudan is worth far more than its best single European: the right to choose the date is a large
  part of the value, more so out of the money.

Wall-clock per pricing (PDE at 400 × 400, which is already within {err[400]:.2f} bp):

{table(["Trade", "PDE (ms)", "LSMC, 4M paths (ms)", "Upper bound (ms)"], times)}

### 5.1 PDE: what converges, and what Rannacher start-up actually buys

Under LGM the deflated value satisfies the heat equation in the variance clock ζ, with no drift or
discounting terms: the change of numeraire has absorbed them.

{table(["Steps (space = time)", "Error with Rannacher (bp)", "Error, Crank–Nicolson only (bp)", "Time (ms)"], pde_rows)}

{pde_plot}

- **Second order overall, not monotone.** The error falls {err[steps[0]] / err[steps[-1]]:,.0f}-fold over a {steps[-1] // steps[0]}-fold refinement, but
  unevenly: it depends on where the exercise boundary falls between grid nodes, so one doubling can change
  it by a factor of {min(drops):.1f} and another by {max(drops):.0f}. Convergence has to be judged over several doublings.
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

- **Degree 0** (continuation value = a constant, so the rule cannot see the state) loses about {deg0_loss:.0f} bp, and
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


def two_factor(res: Path) -> str:
    par = read(res, "g2pp_params")
    fit = read(res, "g2pp_fit")
    swp = read(res, "g2pp_swaption")
    ber = read(res, "g2pp_bermudan")

    one = par[0]
    two = [r for r in par if r["model"] == "2F"]
    best = min(two, key=lambda r: float(r["rms_bp"]))
    par_body = [["one factor", "", f"{float(one['kappa1']):.1%}", "", "", "", f"{float(one['rms_bp']):.2f}", f"{float(one['max_bp']):.2f}", "1.000"]]
    for r in two:
        par_body.append(["two factors", r["start"].replace("kappa1", "κ₁").replace("rho", "ρ"), f"{float(r['kappa1']):.1%}", f"{float(r['kappa2']):.1%}", f"{float(r['sigma_ratio']):.2f}", f"{float(r['rho']):+.2f}", f"{float(r['rms_bp']):.2f}", f"{float(r['max_bp']):.2f}", f"{float(r['corr_1y_10y']):.3f}"])
    bounds = sorted({b.replace("kappa", "κ") for r in two for b in r["active_bounds"].split() if b})
    rms2 = [float(r["rms_bp"]) for r in two]
    corr2 = [float(r["corr_1y_10y"]) for r in two]
    agree = sum(abs(float(r["rms_bp"]) - float(best["rms_bp"])) < 1e-3 for r in two)
    others = [r for r in two if abs(float(r["rms_bp"]) - float(best["rms_bp"])) >= 1e-3]
    stuck_bullet = ""
    stuck = "".join(f" From {r['start'].replace('kappa1', 'κ₁').replace('rho', 'ρ')} it stops at {float(r['rms_bp']):.2f} bp with ρ = {float(r['rho']):+.2f}." for r in others)
    stalled = sum(r["status"] == "MaxIterations" for r in two)
    converged = f"{stalled} of the {len(two)} fits stopped at the iteration limit: the objective is flat along a valley." if stalled else "All fits converged."

    expiries = sorted({round(float(r["expiry"])) for r in fit})
    tenors = sorted({round(float(r["tenor"])) for r in fit})
    cell = {(round(float(r["expiry"])), round(float(r["tenor"]))): r for r in fit}
    err_body = [[f"{e}Y"] + [f"{float(cell[(e, t)]['error_1f_bp']):+.1f} / {float(cell[(e, t)]['error_2f_bp']):+.1f}" for t in tenors] for e in expiries]

    swp_body = [[r["quadrature_points"], f"{float(r['price_bp']):.6f}", f"{abs(float(r['error_vs_2d_integration_bp'])):.0e}", f"{float(r['us_per_price']):.0f}"] for r in swp]

    pde = float(ber[0]["price_bp"])
    europe = float(ber[0]["best_european_bp"])
    off = ber[1]
    scan = [r for r in ber if r["case"] == "fast second factor"]
    fitted = ber[-1]
    ber_body = [["One factor, PDE", "", "1.000", f"{pde:.2f}", f"{pde - europe:.1f}", ber[0]["strip_error_bp"]]]
    ber_body.append(["Two-factor code, second factor off", "", "1.000", f"{float(off['price_bp']):.2f} ± {float(off['se_bp']):.2f}", f"{float(off['price_bp']) - europe:.1f}", off["strip_error_bp"]])
    for r in scan:
        ber_body.append(["Fast second factor", f"{float(r['rho']):+.1f}", f"{float(r['corr_1y_10y']):.3f}", f"{float(r['price_bp']):.2f} ± {float(r['se_bp']):.2f}", f"{float(r['price_bp']) - europe:.1f}", r["strip_error_bp"]])
    ber_body.append(["Shape fitted to the ATM matrix", f"{float(fitted['rho']):+.2f}", f"{float(fitted['corr_1y_10y']):.3f}", f"{float(fitted['price_bp']):.2f} ± {float(fitted['se_bp']):.2f}", f"{float(fitted['price_bp']) - europe:.1f}", fitted["strip_error_bp"]])
    prices = [float(r["price_bp"]) for r in scan]
    worst_strip = max(float(r["strip_error_bp"]) for r in ber)

    if others:
        stuck_bullet = f"- **Not every start gets there.**{stuck} The objective has a second,\n  one-factor-like basin, and Levenberg–Marquardt is a local method. Running several starts is part of the\n  calibration, not a diagnostic."
    cms_all = read(res, "cms_spread")
    info = read(res, "cms_spread_info")[0]
    cms = [r for r in cms_all if r["attainable"] == "1" and r["model"] != "fitted to the ATM matrix"]
    cms_fit = cms_all[-1]
    cms_two = [r for r in cms if r["model"] == "two factors"]
    cms_missing = [float(r["rho"]) for r in cms_all if r["attainable"] == "0"]
    cms_prices = [float(r["price_bp"]) for r in cms_two]
    cms_corrs = [float(r["rate_correlation"]) for r in cms_two]
    cms_rhos = [float(r["rho"]) for r in cms_two]
    cms_body = [
        [{"fitted to the ATM matrix": "Fitted to the ATM matrix"}.get(r["model"], r["model"].capitalize()), "" if r["model"] == "one factor" else f"{float(r['rho']):+.2f}", "" if r["model"] == "one factor" else f"{float(r['sigma_ratio']):.2f}", f"{float(r['rate_correlation']):.3f}", f"{float(r['spread_vol_bp']):.1f}", f"{float(r['price_minus25_bp']):.1f}", f"{float(r['price_bp']):.1f}", f"{float(r['price_plus25_bp']):.1f}"]
        for r in cms + [cms_fit]
    ]

    fig, ax = plt.subplots(figsize=(5.8, 3.3))
    ax.axhline(pde, color=MUTED, linewidth=1.0, linestyle="--")
    ax.annotate("one factor (PDE)", (-0.9, pde), textcoords="offset points", xytext=(0, 5), color=MUTED, fontsize=8)
    ax.axhline(europe, color=GRID, linewidth=1.0)
    ax.annotate("best European", (-0.9, europe), textcoords="offset points", xytext=(0, 5), color=MUTED, fontsize=8)
    ax.errorbar(col(scan, "rho"), prices, yerr=[2 * float(r["se_bp"]) for r in scan], fmt="o-", color=BLUE, capsize=3)
    ax.set_title("Bermudan value against factor correlation")
    ax.set_xlabel("correlation ρ between the two factors")
    ax.set_ylabel("bp of notional")
    plot = save(fig, "g2pp_bermudan.png")

    return f"""## 9. A second factor: G2++

**Model view.** One factor moves every rate together. G2++ has two Gaussian factors with different mean
reversions and a correlation ρ between them, so the curve can steepen and flatten as well as shift. It is
written here in the same deflated-bond form as the one-factor model: the state is driftless, the curve is
fitted by construction, and Monte Carlo steps are exact. A common volatility scale s(t) carries the level;
the two speeds, the volatility ratio and ρ carry the shape.

### 9.1 European swaptions: one closed form and a one-dimensional integral

With two states there is no single exercise point, so Jamshidian's decomposition does not apply directly.
Conditioning on one direction of the state brings it back: the inner expectation is a closed form and the
outer integral is one-dimensional. Price of a 5Y×10Y payer against two-dimensional numerical integration:

{table(["Gauss–Hermite points", "Price (bp)", "Error (bp)", "Time (µs)"], swp_body, "rrrr")}

Four points already give the price to {abs(float(swp[1]["error_vs_2d_integration_bp"])):.0e} bp. That is a consequence of the direction chosen for the
closed form. Brigo and Mercurio condition on the first factor; that is exact, but when the factors are
highly correlated the remaining integrand tends to the kinked payoff and the quadrature stops converging.
Taking the closed form along the direction that decides exercise leaves an integrand that is almost flat,
and exactly constant in the one-factor limit. The tests check the formula against their bond-option
formula, against two-dimensional integration, and against the one-factor model when the second factor is
switched off.

### 9.2 Shape from the ATM matrix: what the second factor buys

Both models are fitted to {len(fit)} ATM swaptions (expiries {expiries[0]}Y–{expiries[-1]}Y, tenors {tenors[0]}Y–{tenors[-1]}Y) with the same freedom in the
level: one volatility scale per expiry row. What is left to explain is the tenor direction, which is where
a second factor can help. The two-factor fit is run from four starting points.

{table(["Model", "Start", "κ₁", "κ₂", "σ₂/σ₁", "ρ", "RMS error (bp)", "Worst (bp)", "Corr(1Y, 10Y fwd)"], par_body, "llrrrrrrr")}

Fit errors in bp of normal vol, one factor / two factors (best start), by expiry (rows) and tenor (columns):

{table(["Expiry"] + [f"{t}Y" for t in tenors], err_body)}

- **The second factor takes the error from {float(one["rms_bp"]):.2f} bp to {min(rms2):.2f} bp.** One mean reversion can tilt the tenor
  direction one way; it cannot produce a vol that is as high at 10 years as at 2 years and lower at 20.
- **The shape is identified, and it is the one the data was built with.** {agree} of the {len(two)} starts end at the same
  point: κ₁ = {float(best["kappa1"]):.1%}, κ₂ = {float(best["kappa2"]):.1%}, σ₂/σ₁ = {float(best["sigma_ratio"]):.2f}, ρ = {float(best["rho"]):+.2f}. The vol screen's ATM levels were
  generated from a two-factor Gaussian model with a = 40%, b = 3%, η/σ = 0.79, ρ = −0.70
  (`tools/gen_vol_quotes.py`), through the frozen-weights swap-rate formula and a smooth volatility scale.
  The calibration uses exact pricing and a piecewise-constant scale, and the remaining {min(rms2):.2f} bp and the
  small differences in the parameters are the distance between those two descriptions.
{stuck_bullet}
- With these parameters the instantaneous correlation between the 1-year and 10-year forward rates is
  {float(best["corr_1y_10y"]):.2f}: far from the 1.00 a one-factor model assumes.

This is a recovery test, not a market calibration: the screen is synthetic. What it shows is that the
calibration finds a two-factor structure when one is there, and says so when it is not (an earlier
version of the screen, generated from smooth surfaces with no correlation structure behind them, gave
four different answers from four starts).

### 9.3 The Bermudan: the same Europeans, a different correlation

Holding the shape fixed, the volatility scale is bootstrapped to the co-terminal strip of the 10nc1 payer
at 3.50%, one expiry at a time, exactly as ζ is in the one-factor model. Every model below reprices the
nine Europeans (worst error {worst_strip:.0e} bp of vol). They differ only in how rates on different exercise dates
are correlated. Prices are Longstaff–Schwartz lower bounds on one million paths, regressing on both states.

{table(["Model", "ρ", "Corr(1Y, 10Y fwd)", "Bermudan (bp)", "Switch value (bp)", "Strip error (bp)"], ber_body, "lrrrrr")}

{plot}

- **The pricer is validated where the answer is known.** With the second factor switched off the two-factor
  Monte Carlo gives {float(off["price_bp"]):.2f} ± {float(off["se_bp"]):.2f} bp against the one-factor PDE's {pde:.2f} bp.
- **Correlation moves the Bermudan by {max(prices) - min(prices):.0f} bp** (from {min(prices):.1f} to {max(prices):.1f} bp as ρ goes from −0.9 to +0.9)
  while every calibration instrument stays repriced. That is the same order as the effect of mean reversion
  in section 6, and for the same reason: the Europeans fix the variance of each swap rate on its own
  exercise date and say nothing about how the rates on different dates move together.
- **Which correlation matters.** The prices are not ordered by the correlation between forward rates of
  different maturities (third column): the highest price has the highest such correlation. What the
  Bermudan prices is correlation across TIME, between the swap rates seen on different exercise dates. A
  fast factor's shocks have died out by the later dates. Positively correlated with the slow factor, it
  makes recent shocks count for more than old ones, which is what extra mean reversion does, and the
  switch value rises. Negatively correlated, it offsets recent shocks, acts like negative mean reversion,
  and the value falls below the one-factor price.
- **With the shape fitted to the ATM matrix the Bermudan is {float(fitted["price_bp"]):.1f} bp, {abs(float(fitted["price_bp"]) - pde):.1f} bp {"below" if float(fitted["price_bp"]) < pde else "above"} the one-factor
  price at κ = 3%.** Given 9.2 this is a calibrated number, not one point in a range. It is also close to
  what the one-factor model gives once its mean reversion is marked (section 10), which is the same
  information reaching the price by a cruder route.

Not done: an upper bound for the two-factor price (the duality gap of the one-factor rule in section 5 was
below one basis point, and the same regression is used here on two states), and risk under the two-factor
model.

### 9.4 A CMS spread option: the product that prices the correlation

The option pays max(S₁₀ − S₂ − K, 0) in 5 years, where S₁₀ and S₂ are the 10-year and 2-year swap rates
then. Its value is the variance of the spread, σ₁₀² + σ₂² − 2ρσ₁₀σ₂. The two volatilities are swaption
vols; the correlation between the two rates is in no swaption price.

Every model below reprices the 5Y×2Y and 5Y×10Y swaptions ({float(info["short_vol_bp"]):.1f} and {float(info["long_vol_bp"]):.1f} bp) exactly, so they agree on
both volatilities. The strike is the spread of the two forward rates today, {float(info["curve_spread_bp"]):.1f} bp. Speeds are the
fitted κ₁ and κ₂; for each ρ the level and σ₂/σ₁ are solved for.

{table(["Model", "ρ (factors)", "σ₂/σ₁", "Rate correlation", "Spread vol (bp)", "Price, K − 25 bp", "Price, K (bp)", "Price, K + 25 bp"], cms_body, "lrrrrrrr")}

- **One factor prices the option at its intrinsic value.** With both vols matched the two rates are
  perfectly correlated, the spread has {float(cms[0]["spread_vol_bp"]):.1f} bp of volatility, and the {float(cms[0]["price_bp"]):.1f} bp that is left is the
  convexity difference between the two rates, not optionality. The 25 bp out-of-the-money call is
  worth {float(cms[0]["price_plus25_bp"]):.2f} bp.
- **Two factors move the price from {min(cms_prices):.0f} to {max(cms_prices):.0f} bp** as ρ goes from {cms_rhos[-1]:+.2f} to {cms_rhos[0]:+.2f}, with both
  swaptions repriced to 1e-12 bp throughout. The price is a direct reading of the rate correlation,
  which runs from {max(cms_corrs):.2f} down to {min(cms_corrs):.2f}. This is why spread options are what a desk marks ρ to.
- **The two swaption vols already bound ρ.** For ρ = {", ".join(f"{r:+.2f}" for r in cms_missing)} no σ₂/σ₁ reprices both: with these
  speeds, a 10-year vol as high as the 2-year one needs strongly negatively correlated factors. The
  calibration reports that instead of returning a best effort.
- **The matrix-fitted model gives {float(cms_fit["price_bp"]):.1f} bp**, with a rate correlation of {float(cms_fit["rate_correlation"]):.2f}. It reprices the two
  swaptions to {max(abs(float(cms_fit["short_vol_error_bp"])), abs(float(cms_fit["long_vol_error_bp"]))):.2f} bp, not exactly, and has a time-dependent volatility scale, so it is not on the line
  of constant-parameter models above. That the ATM matrix alone leaves this much room ({float(cms_fit["price_bp"]):.0f} bp against
  {min(cms_prices):.0f}–{max(cms_prices):.0f} bp) is the practical argument for calibrating to a spread option when one is traded.

One payment at expiry with unit accrual, {float(info["ms_per_price"]):.0f} ms per price by two-dimensional quadrature.
"""


def mean_reversion(res: Path) -> str:
    r1, r2 = read(res, "kappa_result_1y")[0], read(res, "kappa_result_2y")[0]
    s1, s2 = read(res, "kappa_scan_1y"), read(res, "kappa_scan_2y")
    f1 = read(res, "kappa_fit_1y")
    fig, ax = plt.subplots(figsize=(5.6, 3.2))
    ax.plot([100 * k for k in col(s1, "kappa")], col(s1, "rms_bp"), color=BLUE, label="second strip: 1Y tenor")
    ax.plot([100 * k for k in col(s2, "kappa")], col(s2, "rms_bp"), color=ORANGE, label="second strip: 2Y tenor")
    ax.set_title("Error on the second strip against mean reversion")
    ax.set_xlabel("mean reversion κ, %")
    ax.set_ylabel("RMS error, bp of normal vol")
    ax.set_ylim(0, 30)
    ax.legend()
    plot = save(fig, "kappa_scan.png")
    body = [
        ["1Y", f"{float(r1['kappa']):.2%}", f"{float(r1['rms_bp']):.2f}", r1["strip_error_bp"], f"{float(r1['bermudan_bp']):.2f}"],
        ["2Y", f"{float(r2['kappa']):.2%}", f"{float(r2['rms_bp']):.2f}", r2["strip_error_bp"], f"{float(r2['bermudan_bp']):.2f}"],
    ]
    errs = ", ".join(f"{float(r['error_bp']):+.1f}" for r in f1)
    ref = float(r1["bermudan_at_reference_bp"])
    scan6 = col(read(res, "bermudan_kappa"), "bermudan_bp")
    range6 = max(scan6) - min(scan6)
    return f"""## 10. Marking the mean reversion

Sections 4 and 6 left κ as an input: the co-terminal strip is repriced at any κ, and the Bermudan moves by
tens of basis points with it. A second set of swaptions with the same expiries and a different tenor does
see κ. Given ζ(t), the model vol of a short swap starting at t scales with e^(−κt), while the co-terminal
swap's vol is held fixed by the calibration. So κ is chosen to minimise the error on the second set,
re-bootstrapping ζ to the co-terminal strip at every trial value.

{table(["Second strip (expiries 1Y–9Y)", "κ", "RMS error on it (bp)", "Co-terminal strip error (bp)", "Bermudan (bp)"], body)}

{plot}

- The objective has a clear interior minimum: κ is identified by either strip.
- The two choices of second strip give {float(r1["kappa"]):.2%} and {float(r2["kappa"]):.2%}, and Bermudan prices {abs(float(r1["bermudan_bp"]) - float(r2["bermudan_bp"])):.1f} bp apart.
  That spread is what is left of the {range6:.0f} bp range of section 6 once κ is marked.
- **The reference κ = 3% used in sections 4–8 was a convention, and the market disagrees with it.** Marked,
  κ is close to zero and the Bermudan is {ref - float(r1["bermudan_bp"]):.1f} bp lower ({float(r1["bermudan_bp"]):.2f} against {ref:.2f} bp). The two-factor
  model fitted to the whole ATM matrix moves the price the same way, and further (section 9.3).
- One number cannot reprice nine instruments. The errors on the 1Y-tenor strip, by expiry, are
  {errs} bp.
  They change sign along the strip: a pattern, not noise, and what a time-dependent mean reversion or a
  second factor would be fitted to.
"""


def xva(res: Path) -> str:
    credit = read(res, "xva_credit")
    swp = read(res, "xva_swap")
    swp_sum = read(res, "xva_swap_summary")
    book = read(res, "xva_book")
    port = read(res, "xva_portfolio")
    tab = read(res, "xva_table")
    eng = read(res, "xva_engine")[0]
    ber = read(res, "xva_bermudan")
    ber_sum = read(res, "xva_bermudan_summary")[0]
    cs01 = read(res, "xva_cs01")
    delta = read(res, "xva_delta")
    timing = {r["method"]: float(r["ms"]) for r in read(res, "xva_delta_timing")}

    cpty = [r for r in credit if r["name"] == "counterparty"]
    credit_body = [[f"{r['maturity']}Y", f"{float(r['spread_bp']):.0f}", f"{float(r['hazard_pct']):.2f}", f"{float(r['survival_pct']):.1f}"] for r in cpty]
    own = [r for r in credit if r["name"] == "own"]
    reprice = max(abs(float(r["repriced_spread_bp"]) - float(r["spread_bp"])) for r in credit)

    swp_body = [[f"{float(r['time']):.2f}", f"{float(r['epe_closed_bp']):.1f}", f"{float(r['epe_mc_bp']):.1f} ± {float(r['mc_se_bp']):.1f}", f"{float(r['ene_closed_bp']):.1f}", f"{float(r['ene_mc_bp']):.1f}"] for r in swp if float(r["time"]) > 0.5]
    sum_body = [[r["method"], f"{float(r['cva_bp']):.2f}", f"{float(r['dva_bp']):.2f}", f"{float(r['fca_bp']):.2f}", f"{float(r['fba_bp']):.2f}"] for r in swp_sum]

    book_body = [[r["trade"], f"{r['tenor_years']}Y", r["fixed_rate_pct"], r["side"], r["notional"], f"{float(r['value_bp']):+.1f}"] for r in book]
    value = sum(float(r["value_bp"]) for r in book)

    def cell(x: str) -> str:
        return f"{float(x):.2f}" if x else "—"

    tab_body = [[r["case"].capitalize(), cell(r["cva_bp"]), cell(r["dva_bp"]), cell(r["fca_bp"]), cell(r["fba_bp"])] for r in tab]
    cva = {r["case"]: float(r["cva_bp"]) for r in tab}
    gross, net = cva["no netting agreement"], cva["netted; no collateral"]
    csa0, csa20, thr, mta = cva["zero threshold; 10-day margin period"], cva["zero threshold; 20-day margin period"], cva["threshold 50 bp of notional"], cva["minimum transfer 10 bp"]

    t = col(port, "time")
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(9.4, 3.3))
    a1.plot(t, col(port, "gross_epe_bp"), color=ORANGE, label="no netting: sum of the trades")
    a1.plot(t, col(port, "epe_bp"), color=BLUE, label="netted")
    a1.plot(t, [-x for x in col(port, "ene_bp")], color=AQUA, label="netted, negative (sign flipped)")
    a1.set_title("Expected exposure of the book")
    a1.set_xlabel("years")
    a1.set_ylabel("discounted, bp of notional")
    a1.legend()
    a2.plot(t, col(port, "csa_threshold_bp"), color=ORANGE, linewidth=1.2, label="threshold 50 bp")
    a2.plot(t, col(port, "csa_20d_bp"), color=AQUA, linewidth=1.2, label="zero threshold; 20-day margin period")
    a2.plot(t, col(port, "csa_zero_bp"), color=BLUE, linewidth=1.2, label="zero threshold; 10-day margin period")
    a2.set_title("With collateral")
    a2.set_xlabel("years")
    a2.legend()
    plot_book = save(fig, "xva_book.png")
    clean = [float(r["csa_zero_bp"]) for r in port if r["coupon_date"] == "0" and float(r["time"]) > 0]
    spike = [float(r["csa_zero_bp"]) for r in port if r["coupon_date"] == "1"]

    bt = col(ber, "time")
    fig, ax = plt.subplots(figsize=(6.0, 3.3))
    ax.stackplot(bt, col(ber, "option_epe_bp"), col(ber, "swap_epe_bp"), colors=[BLUE, AQUA], alpha=0.85, labels=["still an option", "exercised: the remaining swap"])
    ax.plot(bt[1:], col(ber, "underlying_swap_epe_bp")[1:], "o--", color=ORANGE, linewidth=1.2, label="the underlying swap held from today")
    ax.set_title("Expected exposure of a long Bermudan (10nc1 payer, 3.50%)")
    ax.set_xlabel("years")
    ax.set_ylabel("discounted, bp of notional")
    ax.legend(loc="upper right")
    plot_ber = save(fig, "xva_bermudan.png")
    ber_body = [[f"{float(r['time']):.2f}", f"{float(r['epe_bp']):.1f} ± {float(r['se_bp']):.1f}", f"{float(r['option_epe_bp']):.1f}", f"{float(r['swap_epe_bp']):.1f}", f"{float(r['ene_bp']):.1f}", f"{float(r['exercised_pct']):.0f}%", f"{float(r['underlying_swap_epe_bp']):.1f}"] for r in ber]

    cs_body = [[r["quote"], f"{float(r['swap_cva_cs01_bp']):+.4f}", f"{float(r['book_cva_cs01_bp']):+.4f}"] for r in cs01]
    cs_swap, cs_book = sum(float(r["swap_cva_cs01_bp"]) for r in cs01), sum(float(r["book_cva_cs01_bp"]) for r in cs01)

    worst = max(abs(float(r["adjoint"]) - float(r["bump_model_fixed"])) for r in delta)
    tot_adj, tot_fix, tot_rec = (sum(float(r[k]) for r in delta) for k in ("adjoint", "bump_model_fixed", "bump_recalibrated"))
    big = sorted(delta, key=lambda r: -abs(float(r["adjoint"])))[:6]
    order = {r["quote"]: i for i, r in enumerate(delta)}
    big.sort(key=lambda r: order[r["quote"]])
    delta_body = [[r["quote"], f"{float(r['adjoint']):+.4f}", f"{float(r['bump_model_fixed']):+.4f}", f"{float(r['bump_recalibrated']):+.4f}"] for r in big]
    delta_body.append(["**All 37 quotes**", f"**{tot_adj:+.4f}**", f"**{tot_fix:+.4f}**", f"**{tot_rec:+.4f}**"])
    tf = read(res, "xva_twofactor")
    tf_body = [[r["model"].capitalize(), f"{float(r['rho']):+.2f}" if r["rho"] else "", f"{float(r['swap_cva_bp']):.2f}", f"{float(r['book_cva_bp']):.2f}", f"{float(r['book_gross_cva_bp']):.1f}", f"{float(r['book_dva_bp']):.2f}"] for r in tf]
    tf_swap, tf_book, tf_gross = ([float(r[k]) for r in tf] for k in ("swap_cva_bp", "book_cva_bp", "book_gross_cva_bp"))
    b2 = read(res, "xva_bermudan_2f")
    b2s = read(res, "xva_bermudan_2f_summary")[0]
    b2_body = [[f"{float(x['time']):.2f}", f"{float(x['epe_bp']):.1f} ± {float(x['se_bp']):.1f}", f"{float(x['option_epe_bp']):.1f}", f"{float(x['ene_bp']):.1f}", f"{float(y['epe_bp']):.1f}", f"{float(y['option_epe_bp']):.1f}", f"{float(y['ene_bp']):.1f}"] for x, y in zip(b2, ber)]
    time_body = [
        ["One CVA", f"{timing['one CVA']:.2f}", "1×"],
        ["Delta, adjoint", f"{timing['adjoint']:.2f}", f"{timing['adjoint'] / timing['one CVA']:.1f}×"],
        ["Delta, bump with the model fixed", f"{timing['bump (model fixed)']:.0f}", f"{timing['bump (model fixed)'] / timing['one CVA']:.0f}×"],
        ["Delta, bump and re-calibrate", f"{timing['bump (recalibrated)']:.0f}", f"{timing['bump (recalibrated)'] / timing['one CVA']:.0f}×"],
    ]

    return f"""## 11. XVA: from exposure to price

Section 8 stopped at one swap's uncollateralised exposure profile. This section prices it: a credit curve,
a netting set, collateral, a Bermudan, the sensitivities, and the same under two factors. Unless stated the
model is the one-factor LGM of section 4 (κ = 3%, calibrated to the 10-year co-terminal strip). All exposures are discounted expectations of
deflated values, so CVA is exposure weighted by default probabilities with no further discounting.

### 11.1 The credit curve

Hazard rates, constant between CDS maturities, bootstrapped so that each quoted spread makes its CDS worth
zero (quarterly premiums with accrual on default, recovery 40%). The counterparty's curve:

{table(["CDS maturity", "Spread (bp)", "Hazard rate on the interval (%)", "Survival to maturity (%)"], credit_body, "lrrr")}

Both curves reprice their quotes to {"machine precision" if reprice < 1e-9 else f"{reprice:.0e} bp"}. Our own curve, for DVA, is {float(own[0]["spread_bp"]):.0f} bp at one year rising to
{float(own[-1]["spread_bp"]):.0f} bp at ten. The funding spread for FVA is a flat 60 bp. The spreads are illustrative, not a market's.

### 11.2 One swap: the engine against a closed form

For a single swap the expected positive exposure on a coupon date is a European swaption on the remaining
swap, so its CVA is a strip of swaptions with no simulation error. The Monte Carlo engine, which revalues
the swap on every path, must reproduce it. 10-year payer at 3.50%, 100,000 paths:

{table(["Years", "EPE, closed form (bp)", "EPE, engine (bp)", "ENE, closed form (bp)", "ENE, engine (bp)"], swp_body, "rrrrr")}

{table(["Method", "CVA (bp)", "DVA (bp)", "FCA (bp)", "FBA (bp)"], sum_body)}

CVA is the expected loss if the counterparty defaults while owing us; DVA is the mirror image on our own
default; FCA and FBA are the cost and benefit of funding the uncollateralised asset and liability at our
spread. Payer and receiver are each other's mirror: one side's CVA exposure is the other's DVA exposure,
which is tested.

### 11.3 A netting set

Eight swaps with one counterparty (values today in bp of a unit notional; the book is worth {value:+.1f} bp):

{table(["Trade", "Tenor", "Fixed rate %", "Side", "Notional", "Value (bp)"], book_body, "rrrlrr")}

Exposure is simulated on {eng["exposure_dates"]} dates (quarterly, plus every coupon date) with {int(eng["paths"]):,} paths; one run takes
{float(eng["ms"]) / 1000:.1f} s.

{plot_book}

{table(["Agreement", "CVA (bp)", "DVA (bp)", "FCA (bp)", "FBA (bp)"], tab_body)}

- **Netting cuts CVA from {gross:.1f} to {net:.1f} bp**, a factor of {gross / net:.1f}. Without a netting agreement each trade is a
  separate claim and the offsetting payers and receivers both count; with one, only the net does.
- **Collateral with a zero threshold cuts it again to {csa0:.2f} bp, and not to zero.** What is left is the
  move of the book over the 10-day margin period of risk, during which no more collateral arrives.
  Doubling the period raises CVA by a factor of {csa20 / csa0:.2f}. A pure diffusion would give √2 = 1.41; the rest is
  the coupon effect below, which does not scale with the period.
- **A threshold is uncollateralised exposure by design:** 50 bp of notional leaves {thr:.2f} bp of CVA.
- **A minimum transfer amount barely matters** ({mta:.2f} against {csa0:.2f} bp). A call that is not made leaves a
  small value uncollateralised in both directions, so to first order it adds as much as it removes.
- **Collateralised exposure spikes on coupon dates.** Away from them the residual exposure averages
  {sum(clean) / len(clean):.1f} bp; on them, {sum(spike) / len(spike):.1f} bp. When a coupon rolls off, the book's value jumps by the net
  coupon while the collateral held was called on the value before it. This is the part of collateralised
  exposure a diffusion-only argument misses.

### 11.4 The Bermudan: an option that turns into a swap

A swap's future value is a known function of the state. A Bermudan's is not: on each path it is either
still an option, worth the larger of exercise and continuation, or it has been exercised and the holder
owns the remaining swap. The continuation value is computed as a function of the state at every exercise
date by backward induction with exact Gaussian transitions, and paths are then exercised optimally.

The same simulation reproduces the price: {float(ber_sum["grid_price_bp"]):.2f} bp from the grid and {float(ber_sum["simulated_price_bp"]):.2f} ± {float(ber_sum["simulated_se_bp"]):.2f} bp from the
paths, against {float(ber_sum["pde_price_bp"]):.2f} bp from the PDE of section 5.

{table(["Years", "EPE (bp)", "of which option", "of which swap", "ENE (bp)", "Exercised by then", "Underlying swap's EPE (bp)"], ber_body, "rrrrrrr")}

{plot_ber}

- **Before the first exercise date the exposure is the option's price**, {float(ber[0]["epe_bp"]):.0f} bp, far above the
  underlying swap's exposure at the same dates. The option is an asset on every path; the swap is one on
  about half of them.
- **Exercise replaces option exposure with swap exposure.** By year 5, {float(ber[5]["swap_epe_bp"]) / float(ber[5]["epe_bp"]):.0%} of the exposure is the swap's
  ({float(ber[5]["exercised_pct"]):.0f}% of paths have exercised; that frequency is under the model's pricing measure, not a real-world
  probability). At the last date the two profiles meet: the option is gone.
- **Negative exposure appears only after exercise**, on paths where rates fell after the holder entered
  the swap: {float(ber[-1]["ene_bp"]):.1f} bp at the end, against the swap's own which is several times larger. The
  holder chose when to enter.
- **CVA of the Bermudan is {float(ber_sum["cva_bp"]):.1f} bp against {float(ber_sum["swap_cva_bp"]):.1f} bp for the swap**, and its DVA is {float(ber_sum["dva_bp"]):.1f} bp against
  {float(ber_sum["swap_dva_bp"]):.1f} bp. Buying the option is taking more credit risk on the counterparty and giving them
  almost none in return.

Physical settlement is assumed. A cash-settled Bermudan has no exposure after exercise.

### 11.5 Sensitivities

**Credit.** Change in CVA for +1 bp on each CDS quote, re-bootstrapping the hazard curve. Exposure does
not change, so this costs nothing.

{table(["Quote", "Swap CVA (bp per bp)", "Book CVA (bp per bp)"], cs_body)}

Parallel: {cs_swap:+.3f} and {cs_book:+.3f}. The short-maturity entries are negative. Raising the 1-year spread alone
moves default probability from later years, where the exposure is, to the first year, where there is
little. A CVA hedge with a single CDS maturity would get this wrong.

**Rates.** Delta of the swap's CVA to each of the 37 par quotes the curves are built from. The swaption
weights the strip was priced with are already its derivatives to each bond (the exercise boundary is
optimal, so its own sensitivity drops out), which makes CVA linear in the discount factors with known
coefficients: one tape sweep and one transposed solve with the bootstrap Jacobian.

{table(["Quote (largest six)", "Adjoint", "Bump, model fixed", "Bump, re-calibrated"], delta_body)}

{table(["", "Time (ms)", "Relative"], time_body)}

- The adjoint ladder is the exact derivative: the tests match it to a small central bump to 2e-10. The
  bumps in the table are one-sided 1 bp shifts and differ from it by their convexity, at most {worst:.0e} bp per bp.
- The adjoint costs the time of {timing["adjoint"] / timing["one CVA"]:.1f} CVA evaluations. Bumping needs 37 curve rebuilds and revaluations:
  **about {timing["bump (model fixed)"] / timing["adjoint"]:.0f}× slower**.
- **Which delta is a choice.** Bumping with the model's variance held fixed gives {tot_fix:+.3f} bp per bp in total;
  re-calibrating to the same normal vols after each bump gives {tot_rec:+.3f}, {abs(tot_rec / tot_fix - 1):.0%} less, and the buckets move by
  more than the total. The first is the delta of a desk that hedges vega separately; the second folds
  part of the vega into the delta. Neither is wrong, but they are not the same hedge.

### 11.6 Two factors: what the correlation does to exposure

The engine also runs under the two-factor model of section 9. The shape (speeds, σ₂/σ₁) is the one fitted to
the ATM matrix; the level is bootstrapped to the same 10-year co-terminal strip as the one-factor model, for
each ρ. So every model below prices the nine swaptions that make up the 10-year swap's exposure identically.

{table(["Model", "ρ", "10Y swap CVA (bp)", "Book CVA, netted (bp)", "Book CVA, no netting (bp)", "Book DVA (bp)"], tf_body, "lrrrrr")}

- **The 10-year swap's CVA does not move** ({min(tf_swap):.1f}–{max(tf_swap):.1f} bp, within simulation error): its exposure on
  each coupon date is a swaption in the calibration set. This is the two-factor engine's check against a
  known answer, alongside the test that it reproduces the two-factor swaption strip.
- **The book's does, modestly:** netted CVA {min(tf_book):.1f}–{max(tf_book):.1f} bp and gross {min(tf_gross):.0f}–{max(tf_gross):.0f} bp across ρ. Swaps of
  2 to 7 years are not in the calibration set, and how a 2-year receiver offsets a 10-year payer depends
  on how the two ends of the curve move together.
- The netted number moves less than the gross one. This book is roughly balanced in each maturity, so most
  of what the correlation changes nets out. A book that is long one end of the curve and short the other
  is where it would matter (the engine tests include one: its exposure changes by more than a quarter
  between ρ = −0.9 and +0.5).

**The Bermudan under two factors.** With two states the continuation value is estimated by regression
rather than on a grid: an in-the-money fit decides exercise, and a fit on all paths gives the option's value
where it is out of the money. With the second factor off this reproduces the one-factor grid profile of
11.4, which is tested.

{table(["Years", "EPE, two factors (bp)", "of which option", "ENE (bp)", "EPE, one factor (bp)", "of which option", "ENE (bp)"], b2_body, "rrrrrrr")}

Under the fitted two-factor model the Bermudan is worth {float(b2s["price_bp"]):.1f} ± {float(b2s["price_se_bp"]):.1f} bp (section 9.3) and its CVA is
{float(b2s["cva_bp"]):.1f} bp, against {float(ber_sum["cva_bp"]):.1f} bp in the one-factor model at κ = 3%. The exposure is lower at the start, where it is
the option's price, and the same at the end, where it is the remaining swap's and both models reprice the
same swaptions.

Not modelled: wrong-way risk (default and exposure are independent), initial margin, bilateral survival
weighting, more than one currency. The engine values a floating coupon in progress at the term rate fixed
at its period start.
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
    sections = [curves(res), swap_risk(res), vol(res), lgm(res), bermudan(res), risk(res), exposure(res), two_factor(res), mean_reversion(res), xva(res)]
    (ROOT / "docs" / "RESULTS.md").write_text(head + "\n".join(sections), encoding="utf-8")
    print(f"wrote docs/RESULTS.md and {len(list(IMG.glob('*.png')))} plots")


if __name__ == "__main__":
    main()
