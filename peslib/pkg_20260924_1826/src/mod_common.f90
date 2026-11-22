module mod_common
  implicit none
  ! Centralized common variables to avoid COMMON block layout mismatches
  ! Masses
  double precision, public :: mass(3) = 0.d0
  double precision, public :: mtot = 0.d0
  double precision, public :: mred = 0.d0

  ! Array sizing / channels
  integer, public :: mro = 0
  integer, public :: mvi = 0
  integer, public :: nvi = 0
  integer, public :: n = 0
  integer, public :: nchan = 0
  integer, public :: mro_override = 0
  integer, public :: mvi_override = 0
  integer, public :: nvi_override = 0

  ! Scale factors
  double precision, public :: scales(3) = 0.d0
  double precision, public :: rmlmda = 0.d0

  ! Energies / inputs
  double precision, public :: enrg = 0.d0
  double precision, public :: dnrg = 0.d0
  integer, public :: nnrg = 0
  double precision, public :: emax = 0.d0
  integer, public :: mtr = 0

  ! Ranges
  double precision, public :: rmin = 0.d0
  double precision, public :: rmax = 0.d0
  double precision, public :: smax = 0.d0
  double precision, public :: rdat = 0.d0
  double precision, public :: rmin_override = 0.d0
  double precision, public :: smax_override = 0.d0
  double precision, public :: smin_asym_override = 0.d0

  ! PES-specific max distances (bohr) and thresholding
  double precision, public :: pes_rmax(3) = 0.d0
  double precision, public :: pes_threshold_factor = 2.0d0
  double precision, public :: pes_threshold_w = 2.0d0
  ! Pairwise C6 coefficients for analytic long-range model (order matches rvec/q)
  double precision, public :: pes_C6(3) = 0.d0

  ! Rotors / other integers
  integer, public :: jtot = 0
  integer, public :: ipar = 0
  integer, public :: jpar = 0
  integer, public :: ireact = 1
  integer, public :: jmax = 0
  integer, public :: kmin = 0
  integer, public :: kmax = 0

  ! Output
  integer, public :: nout = 0
  integer, public :: jout = 0

  ! Energy scanning mode (iscale):
  !   0 = linear: enrg, enrg+dnrg, enrg+2*dnrg, ... (default, original ABC)
  !   1 = log-spaced collision energy: ecol_min to ecol_max in nnrg log steps
  !       (collision energy = total energy - threshold)
  integer, public :: iscale = 0
  double precision, public :: ecol_min = 0.d0   ! min collision energy (eV)
  double precision, public :: ecol_max = 0.d0   ! max collision energy (eV)

end module mod_common
