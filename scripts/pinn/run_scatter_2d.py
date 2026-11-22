#!/usr/bin/env python3
"""Compatibility wrapper: the 2D scattering benchmark moved to the converged
CAP reference.  Run scripts/pinn/run_scatter_2d_cap_benchmark.py.

The old no-CAP comparison (this file's former content) compared against a
box reference that is NOT convergent (see src/conventional/scatter_2d.py and
results/convergence_2d.json) and is superseded; it is not kept as a parallel
entry point (single source of truth for the 2D benchmark).
"""
import runpy
import sys
from pathlib import Path

if __name__ == '__main__':
    _target = Path(__file__).resolve().parent / 'run_scatter_2d_cap_benchmark.py'
    sys.argv[0] = str(_target)
    runpy.run_path(str(_target), run_name='__main__')
