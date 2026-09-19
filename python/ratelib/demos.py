"""Locate and run the demo executables (cpp/examples) and parse their CSV (or raw) output."""
from __future__ import annotations

import io
import os
import subprocess
from pathlib import Path

import pandas as pd


def repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def find_build_dir(preferred: str | None = None) -> Path:
    """The CMake build directory. RL_BUILD_DIR env wins; else the first existing preset under build/."""
    if env := os.environ.get("RL_BUILD_DIR"):
        return Path(env)
    build = repo_root() / "build"
    candidates = [preferred] if preferred else []
    candidates += ["dev", "msvc-release", "linux-gcc-release", "msvc-debug", "linux-gcc-san"]
    for c in candidates:
        if c and (build / c).is_dir():
            return build / c
    raise FileNotFoundError(f"no CMake build directory under {build}; configure and build first (see README)")


def demo_path(name: str, build_dir: Path | None = None) -> Path:
    build_dir = build_dir or find_build_dir()
    exe = build_dir / "cpp" / "examples" / (name + (".exe" if os.name == "nt" else ""))
    if not exe.exists():
        raise FileNotFoundError(f"{exe} not built; run `cmake --build --preset <preset>`")
    return exe


def run_demo_text(name: str, *args: object, build_dir: Path | None = None) -> str:
    """Run a demo executable and return its stdout verbatim (for commands that print JSON, e.g. `rl_demo_vol cube`).

    Raises RuntimeError with the C++ error text when the command fails.
    """
    exe = demo_path(name, build_dir)
    env = dict(os.environ)
    env.setdefault("RL_DATA_DIR", str(repo_root() / "data"))
    proc = subprocess.run([str(exe), *map(str, args)], capture_output=True, text=True, env=env)
    if proc.returncode != 0:
        raise RuntimeError(f"{name} failed (exit {proc.returncode}): {proc.stderr.strip()}")
    return proc.stdout


def run_demo(name: str, *args: object, build_dir: Path | None = None) -> pd.DataFrame:
    """Run a demo executable and parse its CSV stdout into a DataFrame.

    Raises RuntimeError with the C++ error text (e.g. 'not implemented: Date::daysFromCivil') when the
    stubs are still unimplemented — that is the notebook's signal that there is work left to do.
    """
    return pd.read_csv(io.StringIO(run_demo_text(name, *args, build_dir=build_dir)))
