MODULE OneDim_Parms_Module
   USE Numeric_Kinds_Module
   IMPLICIT NONE
   SAVE
   INTEGER  IArran
   INTEGER, PARAMETER :: nmax=2**10  ! nmax+1=number of radial grid points (increase to get convergence)
   INTEGER, PARAMETER :: MaxYterms=5! Minimum yterms (usually 3)
   INTEGER, PARAMETER :: ndaf=2     ! Order of differential equation (always 2)
   INTEGER, PARAMETER :: mval=32    ! order of Hermite DAF expansion (always 32)
   
   LOGICAL FitDAFEnergies           ! IF TRUE the DAF eigenenergies will be fit
   LOGICAL FitNumerovEnergies       ! IF TRUE the Numerov eigenenergies will be fit
   LOGICAL FitFourierGridEnergies   ! IF TRUE the FourierGrid eigenenergies will be fit
   LOGICAL Continuum                ! IF TRUE Scattering eigenfunctions will be calculated
   INTEGER jmin                     ! Minimum rotational state (usually 0)
   INTEGER jmax                     ! Maximum rotational state
   INTEGER MinElec                  ! Minimum electronic state (usually 0)
   INTEGER MaxElec                  ! Maximum electronic state (usually 0 or 1)
   INTEGER NEnergy                  ! Number of continuum energies to calculate
   REAL(KIND=WP_Kind)  Emin         ! Minimum continuum energy (Hartree Units)
   REAL(KIND=WP_Kind)  Emax         ! Maximum continuum energy (Hartree Units)
   REAL(KIND=WP_Kind)  MassA        ! Atomic mass of atom A     
   REAL(KIND=WP_Kind)  MassB        ! Atomic mass of atom A
   REAL(KIND=WP_Kind)  ReducedMass  ! System reduced mass
END MODULE OneDim_Parms_Module