SUBROUTINE Constant
!
!-----------------------------------------------------------------
!  This routine supplies constants and conversion factors.
!-----------------------------------------------------------------
!  pi            3.14159265358979323846264438327950288......
!  twopi         2*pi
!  halfpi        pi/2
!  fourthpi      pi/4
!-----------------------------------------------------------------
!  Conversion factors
!-----------------------------------------------------------------
!  radtodeg      converts radians to degrees
!  autoev        converts hartrees to electron volts.
!  amutoau       atomic mass units to atomic units.
!  angtoau       angstroms to bohr.
!  cmm1toau     inverse centimeters to atomic units
!-----------------------------------------------------------------
!  Physical constants (1986 values, Phys. Today, Aug. 1987)
!-----------------------------------------------------------------
!  avagadro = 6.0221367d+23     molecules/mole
!  clight   = 2.99792458d+10    cm/sec (defined constant)
!  boltzman = 1.380658d-16      ergs/Kelvin
!-----------------------------------------------------------------
!  prefixes
!-----------------------------------------------------------------
!  exa      = 1.d+18
!  peta     = 1.d+15
!  tera     = 1.d+12
!  giga     = 1.d+9
!  mega     = 1.d+6
!  kilo     = 1.d+3
!  hecto    = 1.d+2
!  deca     = 1.d+1
!  deci     = 1.d-1
!  centi    = 1.d-2
!  milli    = 1.d-3
!  micro    = 1.d-6
!  nano     = 1.d-9
!  pico     = 1.d-12
!  femto    = 1.d-15
!  atto     = 1.d-18
!-----------------------------------------------------------------
!  fuzz is a error factor for comparing real numbers
!-----------------------------------------------------------------
USE Numeric_Kinds_Module
USE FileUnits_Module
USE Fuzzy_Module
USE Convrsns_Module
USE Constnts_Module
USE Prefix_Module
USE Time_Conv_Module
USE Numbers_Module
IMPLICIT NONE
SAVE
!-----------------------------------------------------------------
!  Conversion factors
!-----------------------------------------------------------------
radtodeg = fortyfive/fourthpi
amutoau  = 1822.888506d0
angtoau  = 0.529177249d0
cmm1toau = 4.5563352672d-6
!-----------------------------------------------------------------
!  Physical  constants
!-----------------------------------------------------------------
avagadro = 6.0221367d+23
clight   = 2.99792458d+10
boltzman = 1.380658d-16
hartree = 4.3597482d-18
hbar    = 1.05457266d-34
!-----------------------------------------------------------------
!  fuzz
!-----------------------------------------------------------------
fuzz = 1.0d-10
!
! Time Conversion. Converts from Atomic units to femto-seconds.
!
time_conv = hbar/hartree/femto

RETURN
ENDSUBROUTINE Constant
