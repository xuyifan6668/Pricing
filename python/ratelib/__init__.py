"""ratelib — Python layer of the rates pricing library.

Helpers for validating market data and for driving the C++ executables from notebooks. There are no
Python bindings to the C++ library itself (docs/DESIGN.md).
"""
from __future__ import annotations

from ratelib.demos import find_build_dir, run_analysis, run_demo, run_demo_text

__all__ = ["find_build_dir", "run_analysis", "run_demo", "run_demo_text"]
__version__ = "0.1.0"
