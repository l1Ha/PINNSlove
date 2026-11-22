"""Physics-Informed Neural Network (PINN) Quantum Solvers Subpackage.

Authentic Neural Solvers:
------------------------
- VariationalPINN: Deep-Ritz / Sobolev-variational physics-informed neural network
  for bound states in molecular potentials (1D Morse, 2D/3D Water).
- ScatteringPINN1D: Physics-informed neural network for 1D fixed-energy scattering.
"""
from .pinn_bound import VariationalPINN
from .pinn_scatter import ScatteringPINN1D
from .pinn_scatter_2d import ScatteringPINN2D

__all__ = [
    'VariationalPINN',
    'ScatteringPINN1D',
    'ScatteringPINN2D',
]
