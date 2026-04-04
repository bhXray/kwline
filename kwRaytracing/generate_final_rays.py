#!/usr/bin/env python3
"""Generate only the ray files used by plot.ipynb."""
import os
import shutil
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parent
EXE = ROOT / "trace_one_ray"
OUT_DIR = ROOT / "output"
TARGET_DIR = ROOT / "scan_rays"

# (alpha, beta, target filename)
TARGETS = [
    (-18.00, -5.0, "ray_a-18.00_b-5.0.txt"),
    (-15.00, -5.0, "ray_a-15.00_b-5.0.txt"),
    (-14.00, -5.0, "ray_a-14.00_b-5.0.txt"),
    (-13.1407, -5.0, "ray_a-13.1407_b-5.0.txt"),
    (-13.0402, -5.0, "ray_a-13.0402_b-5.0.txt"),
    (-13.7437, -5.0, "ray_a-13.7437_b-5.0.txt"),
    (-10.10, -7.0, "ray_a-10.10_b-7.0.txt"),
    (-12.7387, -7.0, "ray_a-12.7387_b-7.0.txt"),
    (-11.70, -10.0, "ray_a-11.70_b-10.0.txt"),
]

# Fixed physics parameters used throughout this project.
A_SPIN = 0.998
LAMBDA_BAR = 0.90
COBS = 80.0
ROBS = 200.0
SCAL = 1.0


def run_one(alpha: float, beta: float, target_name: str) -> None:
    line = f"{alpha:.8f} {beta:.6f} {A_SPIN:.6f} {LAMBDA_BAR:.6f} {COBS:.6f} {ROBS:.6f} {SCAL:.6f}\n"
    subprocess.run(
        [str(EXE)],
        input=line,
        text=True,
        cwd=str(ROOT),
        check=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )

    ray1 = OUT_DIR / "ray1.txt"
    if not ray1.exists():
        raise FileNotFoundError("trace_one_ray did not produce output/ray1.txt")

    TARGET_DIR.mkdir(exist_ok=True)
    shutil.move(str(ray1), str(TARGET_DIR / target_name))


if __name__ == "__main__":
    if not EXE.exists():
        raise FileNotFoundError("trace_one_ray not found. Run 'make one' first.")

    for alpha, beta, name in TARGETS:
        run_one(alpha, beta, name)
        print(f"saved scan_rays/{name}")
