SUBROUTINE Time_Delays ( NOpen, S_Real, S_Imag, SE_real, SE_imag, TD_Real, TD_Imag)
USE Numeric_Kinds_Module
USE FileUnits_Asymptotic_Module
USE physical_constants_Module
IMPLICIT NONE
CHARACTER(LEN=11), PARAMETER:: ProcName='Time_Delays'
CHARACTER(LEN=6) Print_Flag
!
!     Author: Gregory A. Parker, Department of Physics and Astronomy, University of Oklahoma
!
!     This SUBROUTINE calculates the time delay defined as:
!     deltime(IState,JState) = Re[-ihbar(1.0d0/S(IState,JState))*dS(IState,JState)/dE]
!     where S is the S-Matrix and the derivative is with respect to the total energy (E)
!
!     References
!     Z. Darakjian, E.F. Hayes, G.A. Parker, E.A. Butcher and J.D. Kress JCP 95, 2516-2522 (1991)
!     F.T. Smith, Phys. Rev. 118, 349 (1960) [Original definition of Time-Delays]
!
!     Required Input <=====
!     NOpen           Size of ALL matrices (NOpen by NOpen)
!     S_Real      Real part of the S-Matrix
!     S_Imag      Imaginary part of the S-Matrix
!     SE_Real     Energy derivative of the real part of the S-Matrix
!     SE_Imag     Energy derivative of the imaginary part of the S-Matrix
!                 The units of Energy must be Hartree Atomic Units.
!
!     On return  =====>
!     TD_Real     Real part of the Time-Delays (TD)
!     TD_Imag     Imaginary part of Time-Delays (TD)
!                 Units of the lifetime matrix are Femto-seconds
!
!     This routine is called by:
!         Asymptotic
!     This routine calls:
!         Matrix_Out
!
INTEGER, INTENT(IN):: NOpen
REAL(Kind=WP_Kind), INTENT(IN)::   S_Real(NOpen,NOpen),  S_imag(NOpen,NOpen)  ! S-Matrix
REAL(Kind=WP_Kind), INTENT(IN)::  SE_Real(NOpen,NOpen), SE_Imag(NOpen,NOpen)  ! Energy derivative of S-Matrix
REAL(Kind=WP_Kind), INTENT(OUT):: TD_Real(NOpen,NOpen), TD_Imag(NOpen,NOpen)  ! Time-Delays

WRITE(Out_Unit,*)
WRITE(Out_Unit,*)
WRITE(Out_Unit,*)
WRITE(Out_Unit,*)'Entering:', ProcName
CALL PoptAsy(ProcName, Print_Flag)

TD_Real=mhbar*(S_Real*SE_Imag-S_Imag*SE_Real)/(S_Real**2 + S_Imag**2)     ! Calculate real part of time-delays
TD_Imag=mhbar*(S_Real*SE_Imag+S_Imag*SE_Real)/(S_Real**2 + S_Imag**2)     ! Calculate imaginary part of time-delays

CALL Matrix_Out(TD_Real,NOpen,NOpen,'TD_Real', 'Real Part of the Time-delays', Print_Flag)       ! Print real part
CALL Matrix_Out(TD_Imag,NOpen,NOpen,'TD_Imag', 'Imaginary Part of the Time-delays', Print_Flag)  ! Print imaginary part

WRITE(Out_Unit,*)'Leaving:', ProcName
RETURN
ENDSUBROUTINE Time_Delays
