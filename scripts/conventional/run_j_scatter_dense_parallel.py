#!/usr/bin/env python3
"""Dense 121-point parallel runner for J > 0 coupled scattering on BKMP2.

Evaluates J in [0, 1, 2, 3], parities [+1, -1], and 121 collision energies
across all CPU cores in parallel, generating high-density Figure 12 (> 100 points per curve).
"""
from __future__ import annotations

import json
import multiprocessing as mp
import os
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.j_scatter_comprehensive import (
    solve_one_point, plot_all_results, ENERGIES, JS
)


def task_worker(args):
    J, p, Ec = args
    try:
        return solve_one_point(J=J, parity=p, Ec=Ec, nR=350)
    except Exception as exc:
        print(f"Error J={J}, p={p}, Ec={Ec}: {exc}")
        return None


def main():
    tasks = []
    for J in JS:
        for p in [+1, -1]:
            for Ec in ENERGIES:
                tasks.append((J, p, Ec))

    n_cpus = min(os.cpu_count() or 4, 96)
    print("=" * 80)
    print(f"J > 0 DENSE COUPLED SCATTERING: {len(tasks)} TOTAL TASKS (121 ENERGIES > 100)")
    print(f"Parallelizing across {n_cpus} cores...")
    print("=" * 80)

    t0 = time.time()
    with mp.Pool(processes=n_cpus) as pool:
        results = pool.map(task_worker, tasks)

    results = [r for r in results if r is not None]
    dt = time.time() - t0
    print(f"Finished {len(results)}/{len(tasks)} solves in {dt:.1f} s ({dt/len(results)*1000:.1f} ms/solve).")

    print("\nGenerating comprehensive 9-panel Figure 12 (> 100 points per curve)...")
    plot_all_results(results)
    print("[✓] Saved dense Figure 12: article/fig/fig_j_scatter_comprehensive.png")


if __name__ == '__main__':
    main()
