"""Scripts namespace package for PINNSlove.

Subpackages:
------------
- scripts.pinn: Execution scripts for neural network solvers (1D to 12D).
- scripts.conventional: Execution scripts for classical numerical solvers (DVR, FD, Close-Coupling).
- scripts.validation: Verification, environment checks, and sub-per-mille audit suites.
- scripts.plotting: Publication figure generation and visualization scripts.
"""
import os
import sys
import importlib

_parent = os.path.dirname(__file__)
__path__ = [
    _parent,
    os.path.join(_parent, 'pinn'),
    os.path.join(_parent, 'conventional'),
    os.path.join(_parent, 'validation'),
    os.path.join(_parent, 'plotting'),
]

_SCRIPT_CACHE = {}

def __getattr__(name):
    if name in _SCRIPT_CACHE:
        return _SCRIPT_CACHE[name]
    for sub in ('pinn', 'conventional', 'validation', 'plotting'):
        target = f'{__name__}.{sub}.{name}'
        try:
            mod = importlib.import_module(target)
            _SCRIPT_CACHE[name] = mod
            globals()[name] = mod
            sys.modules[f'{__name__}.{name}'] = mod
            return mod
        except (ModuleNotFoundError, AttributeError):
            pass
    raise AttributeError(f"module '{__name__}' has no attribute {name!r}")
