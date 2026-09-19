"""ratelib — Python layer of the exotic rates pricing library.

ships only helpers for driving the C++ demo executables from notebooks. The pybind11 module
`ratelib._ratelib` arrives  (docs/DESIGN.md).
"""
from __future__ import annotations

from ratelib.demos import find_build_dir, run_demo, run_demo_text

__all__ = ["find_build_dir", "run_demo", "run_demo_text"]
__version__ = "0.1.0"
