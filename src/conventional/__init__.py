"""Conventional (Classical Numerical / Traditional) Quantum Solvers Subpackage.

Provides finite-difference, DVR, close-coupling, Numerov, Lanczos, and
analytic reference solvers for comparison with neural methods.
"""
from .reference import (
    fd_bound_1d,
    scatter_fd_1d,
    LegendreAngularDVR,
    dvr_2d,
    dvr_3d,
    gauss_legendre_1d,
)
from .bkmp2 import (
    V_h3_angstrom,
    asymptotic_energy,
    get_lib,
    build_library,
    BKMP2Unavailable,
)
from .bkmp2_j_coupled import (
    solve_j_scattering_bkmp2,
    compute_bkmp2_potential_matrix,
)
from .diatom_diatom_4d import DiatomDiatom4DScattering
from .diatom_diatom_6d import DiatomDiatom6DScattering
from .hamiltonians import (
    MorseProblem1D,
    Water2DProblem,
    Water3DProblem,
)
from .hamiltonians_nd import CoupledAnharmonicND
from .j_frames import JAngularCoupling
from .scatter_2d import (
    assemble_2d_scatter,
    scatter_fd_2d,
)
from .scatter_extract import (
    fd2_4th_dirichlet,
    fit_wavevector,
    extract_smatrix,
    solve_radial_single,
    free_incoming,
    free_outgoing,
)
from .semiclassical import SemiclassicalLanger1D

__all__ = [
    'fd_bound_1d',
    'scatter_fd_1d',
    'LegendreAngularDVR',
    'dvr_2d',
    'dvr_3d',
    'gauss_legendre_1d',
    'V_h3_angstrom',
    'asymptotic_energy',
    'get_lib',
    'build_library',
    'BKMP2Unavailable',
    'solve_j_scattering_bkmp2',
    'compute_bkmp2_potential_matrix',
    'DiatomDiatom4DScattering',
    'DiatomDiatom6DScattering',
    'MorseProblem1D',
    'Water2DProblem',
    'Water3DProblem',
    'CoupledAnharmonicND',
    'JAngularCoupling',
    'assemble_2d_scatter',
    'scatter_fd_2d',
    'fd2_4th_dirichlet',
    'fit_wavevector',
    'extract_smatrix',
    'solve_radial_single',
    'free_incoming',
    'free_outgoing',
    'SemiclassicalLanger1D',
]
