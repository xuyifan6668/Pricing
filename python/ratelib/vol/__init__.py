"""Swaption vol cube: screen validation, cube assembly orchestration and the calibration report.

docs/DESIGN.md §3: the per-node SABR fit kernel is C++ (called ~1e4 times per cube-vega ladder); the loop
over nodes, the beta/shift policy, validation and the report are Python because they change weekly.
`vol_quotes` checks the screen before any fit; `cube_builder` drives `rl_demo_vol`; `report` renders the
residual table and the negative-density map.
"""
from ratelib.vol.cube_builder import USD_POLICY, CubeBuild, SabrPolicy, arbitrage_scan, build, build_cube, fit_nodes, fit_residuals, query_smile
from ratelib.vol.report import CubeReport, load_arbitrage_csv, load_fit_csv
from ratelib.vol.vol_quotes import VolNode, VolQuotes, VolQuoteValidationError, load_vol_quotes, tenor_years, validate_vol_quotes

__all__ = [
    "VolNode",
    "VolQuotes",
    "VolQuoteValidationError",
    "load_vol_quotes",
    "validate_vol_quotes",
    "tenor_years",
    "SabrPolicy",
    "USD_POLICY",
    "CubeBuild",
    "fit_nodes",
    "fit_residuals",
    "arbitrage_scan",
    "build_cube",
    "query_smile",
    "build",
    "CubeReport",
    "load_fit_csv",
    "load_arbitrage_csv",
]
