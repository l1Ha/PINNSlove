"""Backwards compatibility shim for src.pinn.pinn_scatter_2d."""
import sys
from src.pinn.pinn_scatter_2d import *  # noqa: F401,F403
from src.pinn import pinn_scatter_2d as _target
for _k, _v in _target.__dict__.items():
    if not _k.startswith("__"):
        globals()[_k] = _v
