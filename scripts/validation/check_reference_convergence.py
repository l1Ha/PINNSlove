"""Grid-convergence check of the DVR reference solvers (2D and 3D)."""
import sys, time
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from src.potentials import WaterPES
from src.reference import dvr_2d, dvr_3d

pes = WaterPES()
print('2D DVR grid convergence (E - Vmin):', flush=True)
prev = None
for (nR, ng) in [(300, 48), (500, 80), (800, 120)]:
    t0 = time.time()
    v, _, _, _, _, _ = dvr_2d(pes, pes.re, 0.55, 3.2, nR, ng, k=5)
    d = '' if prev is None else f'  max|diff| vs prev = {np.max(np.abs(v-prev)):.2f}'
    print(f'  nR={nR:4d} n_gamma={ng:3d}: {np.round(v-pes.V_min(), 3)}{d}  ({time.time()-t0:.0f}s)', flush=True)
    prev = v

print('3D DVR grid convergence (E - Vmin):', flush=True)
prev = None
for (nr, nR, ng) in [(60, 66, 30), (90, 100, 44), (120, 132, 56)]:
    t0 = time.time()
    v, _, _, _ = dvr_3d(pes, 0.55, 2.4, nr, 0.55, 3.0, nR, ng, k=4)
    d = '' if prev is None else f'  max|diff| vs prev = {np.max(np.abs(v-prev)):.2f}'
    print(f'  nr={nr:3d} nR={nR:3d} ng={ng:2d}: {np.round(v-pes.V_min(), 3)}{d}  ({time.time()-t0:.0f}s)', flush=True)
    prev = v
