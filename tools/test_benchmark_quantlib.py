"""The QuantLib comparison as a test. Skipped where QuantLib or the C++ build is not available.

    uv run --with QuantLib --with pytest pytest tools/test_benchmark_quantlib.py
"""
from __future__ import annotations

import pytest

ql = pytest.importorskip("QuantLib", reason="QuantLib is needed only for the independent benchmark")

import benchmark_quantlib  # noqa: E402
from ratelib.demos import find_build_dir, run_analysis  # noqa: E402


@pytest.fixture(scope="module")
def rows():
    try:
        build = find_build_dir()
        run_analysis("benchmark", build_dir=build)
    except FileNotFoundError as e:
        pytest.skip(f"no C++ build to benchmark: {e}")
    # Coarser lattice and grid than the published table: the test bounds the disagreement, it does not measure it.
    return benchmark_quantlib.compare(build / "notebook-results" / "benchmark", tree_steps=400, fd_time=300, fd_space=300)


def worst(rows, prefix):
    found = [abs(r["difference"]) for r in rows if r["quantity"].startswith(prefix)]
    assert found, f"no rows for {prefix}"
    return max(found)


def test_every_case_and_quantity_is_compared(rows):
    assert len({r["case"] for r in rows}) == 4
    assert sum(r["quantity"].startswith("European") for r in rows) == 4 * 9


def test_swap_value_and_par_rate_agree_to_rounding(rows):
    # Same discount factors on the same dates: only accruals and the leg arithmetic are independent.
    assert worst(rows, "swap value") < 1e-6
    assert worst(rows, "par rate") < 1e-6


def test_european_swaptions_agree_with_the_hull_white_closed_form(rows):
    # LGM with Jamshidian's decomposition here; Hull-White bond options summed in QuantLib.
    assert worst(rows, "European") < 1e-3


def test_bermudan_agrees_with_two_other_numerical_methods(rows):
    # Crank-Nicolson in the variance clock here; a trinomial lattice and a finite-difference solver in the short rate there.
    assert worst(rows, "Bermudan, PDE against QuantLib lattice") < 0.5
    assert worst(rows, "Bermudan, PDE against QuantLib finite differences") < 0.2
