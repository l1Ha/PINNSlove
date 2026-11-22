"""Version and environment compatibility shims.

Handles:
* NumPy < 2.0 (no ``np.trapezoid``; falls back to ``np.trapz``)
* PyTorch not installed (the reference solvers and the coupled-channel
  machinery are pure NumPy/SciPy and must remain importable)
* SciPy < 1.10 (``scipy.sparse.block_diag`` ``format=`` keyword)
"""
import numpy as np

# --- NumPy < 2.0 compatibility -------------------------------------------
if hasattr(np, 'trapezoid'):
    trapezoid = np.trapezoid
else:                                   # NumPy 1.x
    trapezoid = np.trapz


# --- optional PyTorch -----------------------------------------------------
def has_torch():
    """True when PyTorch can be imported (needed only for the PINN solvers)."""
    try:
        import torch  # noqa: F401
        return True
    except Exception:
        return False


def get_torch():
    """Return the torch module or raise a clear, actionable error."""
    try:
        import torch
        return torch
    except Exception as exc:            # pragma: no cover
        raise ImportError(
            "PyTorch is required for the PINN solvers but is not available.\n"
            "  pip install torch --index-url https://download.pytorch.org/whl/cpu\n"
            "(The classical reference solvers in src.reference, src.j_frames "
            "and src.scatter_extract work without PyTorch.)"
        ) from exc


def describe_environment():
    """Return a dict with the versions of the libraries that matter."""
    import platform
    import scipy
    info = {
        'python': platform.python_version(),
        'numpy': np.__version__,
        'scipy': scipy.__version__,
        'torch': None,
        'cuda': False,
    }
    try:
        import torch
        info['torch'] = torch.__version__
        info['cuda'] = bool(torch.cuda.is_available())
    except Exception:
        pass
    try:
        import matplotlib
        info['matplotlib'] = matplotlib.__version__
    except Exception:
        info['matplotlib'] = None
    return info
