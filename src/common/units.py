"""Unit system: length [Angstrom], mass [amu], energy [cm^-1], angle [radian]."""
import numpy as np
from scipy import constants as const

# 1 Joule expressed in cm^-1
J_TO_CM1 = 1.0 / (const.h * const.c * 100.0)

# hbar^2 / (2 * 1amu * (1 Angstrom)^2)  in cm^-1 * A^2  (~16.8576)
H2_OVER_2AMU = (const.hbar ** 2 / (2.0 * const.atomic_mass * 1e-20)) * J_TO_CM1


def kinetic_factor(mass_amu):
    """hbar^2/(2m) in cm^-1*A^2 for a mass given in amu (T = -B d^2/dx^2)."""
    return H2_OVER_2AMU / np.asarray(mass_amu, dtype=float)


# Atomic masses (most abundant isotopes), amu
M_O = 15.994915
M_H = 1.007825
