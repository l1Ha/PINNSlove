"""PINNSlove: Quantum Dynamics with Physics-Informed Neural Networks and Conventional Methods.

Directory Architecture:
-----------------------
- src.pinn: Physics-Informed Neural Network solvers for molecular quantum systems.
- src.conventional: Classical numerical reference solvers (DVR, FD, Close-Coupling, Numerov).
- src.common: Physical units, compatibility shims, and shared potentials.
"""
import os
import sys
import importlib

_parent = os.path.dirname(__file__)
__path__ = [
    _parent,
    os.path.join(_parent, 'pinn'),
    os.path.join(_parent, 'conventional'),
    os.path.join(_parent, 'common'),
]

# Explicit subpackages
from . import pinn
from . import conventional
from . import common

# Convenience direct exports from subpackages
from .pinn import *
from .conventional import *
from .common import *

# Module-level __getattr__ for backwards-compatible module imports (e.g. from src import bkmp2)
_SUB_MODULE_CACHE = {}

def __getattr__(name):
    if name in _SUB_MODULE_CACHE:
        return _SUB_MODULE_CACHE[name]
    for sub in ('pinn', 'conventional', 'common'):
        target = f'{__name__}.{sub}.{name}'
        try:
            mod = importlib.import_module(target)
            _SUB_MODULE_CACHE[name] = mod
            globals()[name] = mod
            sys.modules[f'{__name__}.{name}'] = mod
            return mod
        except (ModuleNotFoundError, AttributeError):
            pass
    raise AttributeError(f"module '{__name__}' has no attribute {name!r}")
