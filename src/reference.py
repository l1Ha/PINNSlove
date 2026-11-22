"""Backwards compatibility shim for src.conventional.reference."""
import sys
from src.conventional import reference as _target

# Re-export everything including private symbols
for _k, _v in _target.__dict__.items():
    if not _k.startswith('__'):
        globals()[_k] = _v

# Also copy __all__ if present
if hasattr(_target, '__all__'):
    __all__ = _target.__all__
