"""Notebook helper: make the in-repo `ratelib` package importable without installing it."""
import sys
from pathlib import Path

_ROOT = Path(__file__).resolve().parents[1]
for p in (_ROOT / "python", _ROOT / "tools"):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))

from ratelib.demos import find_build_dir, run_analysis, run_demo  # noqa: E402,F401
