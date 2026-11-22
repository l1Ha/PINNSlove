"""One-off upsized-PINN trial for E_c=200 against the converged CAP reference.

Not part of the regression suite: exploratory capacity study.  The final
benchmark remains scripts/pinn/run_scatter_2d_cap_benchmark.py.
"""
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import torch

from src.potentials import WaterPES
from src.pinn.pinn_scatter_2d import ScatteringPINN2D
from run_scatter_2d_cap_benchmark import compare, fd

pes = WaterPES()
E_c = 200.0
chi_ref, R, g, K, E_tot = fd(pes, E_c, 400, 60, 100)
sp = ScatteringPINN2D(pes, E_c, R_hi=12.0, n_quad=(400, 80), hidden=(4, 128),
                      n_fourier=64, seed=0, device='cuda',
                      cap_start=5.0, cap_eta=100.0, cap_width=2.0)
t0 = time.time()
sp.train(steps=10000, lr0=3e-3, verbose_every=2500, lbfgs=True)
ps = sp.evaluate(R, g)
print(f'upsized PINN E_c=200: L2(R<4.9) = '
      f'{compare(ps, R, g, chi_ref, R, g, 4.9):.4f}  L2(R<6.9) = '
      f'{compare(ps, R, g, chi_ref, R, g, 6.9):.4f}  ({time.time()-t0:.0f}s)',
      flush=True)
