SUBROUTINE T_Matrix(NOpen, S_Real, S_Imag, SE_Real, SE_Imag, T_Real, T_Imag, TE_Real, TE_Imag, EDeriv)
USE Numeric_Kinds_Module
USE FileUnits_Asymptotic_Module
IMPLICIT NONE
CHARACTER(LEN=8), PARAMETER:: ProcName='T_Matrix'
CHARACTER(LEN=6) Print_Flag
!
!     Author: Gregory A. Parker, Department of Physics and Astronomy, University of Oklahoma
!
!     This SUBROUTINE calculates the T-Matirx (T) defined by:
!           T = I - S
!
!     References
!     R. T Pack and G.A. Parker, JCP 87, 3888-3921 (1987)
!
!     Required Input <=====
!     NOpen           Size of ALL matrices (NOpen by NOpen)
!     EDeriv=True If ernergy derivatives of the S-Matrix exist
!     S_Real      Real part of the S-Matrix
!     S_Imag      Imaginary part of the S-Matrix
!     SE_Real     Energy derivative of the real part of the S-Matrix
!     SE_Imag     Energy derivative of the imaginary part of the S-Matrix
!
!     On return  =====>
!     T_Real     Real part of the T-Matrix (T)
!     T_Imag     Imaginary part of T-Matrix (T)
!     TE_Real    Energy derivative of the real part of the T-Matrix (T)
!     TE_Imag    Energy derivative of the imaginary part of the T-Matrix (T)
!
!     This routine is called by:
!         Asymptotic
!     This routine calls:
!         Unit_Matrix
!         Matrix_Out
!         Popt
INTEGER, INTENT(IN):: NOpen
LOGICAL, INTENT(IN):: EDeriv
REAL(Kind=WP_Kind), INTENT(IN)::  S_Real(NOpen,NOpen),  S_Imag(NOpen,NOpen)
REAL(Kind=WP_Kind), INTENT(IN)::  SE_Real(NOpen,NOpen), SE_Imag(NOpen,NOpen)
REAL(Kind=WP_Kind), INTENT(OUT):: T_Real(NOpen,NOpen),  T_Imag(NOpen,NOpen)
REAL(Kind=WP_Kind), INTENT(OUT):: TE_Real(NOpen,NOpen), TE_Imag(NOpen,NOpen)
REAL(Kind=WP_Kind), ALLOCATABLE:: Unit_Mat(:,:)

WRITE(Out_Unit,*)
WRITE(Out_Unit,*)
WRITE(Out_Unit,*)
WRITE(Out_Unit,*)'Entering:', ProcName
ALLOCATE(Unit_Mat(NOpen,NOpen))

CALL PoptAsy(ProcName, Print_Flag)

CALL Unit_Matrix(NOpen,Unit_Mat)

!
! Calculate the T_matrix
T_Real=Unit_Mat-S_Real                                                                 ! Calculate Real(T)
CALL Matrix_Out(T_Real,NOpen,NOpen,'T_Real','Real Part of T-Matrix', Print_Flag)       ! Print Real(T)
T_Imag=-S_Imag                                                                         ! Calculate Imag(T)
CALL Matrix_Out(T_Imag,NOpen,NOpen,'T_Imag','Imaginary Part of T-Matrix', Print_Flag)  ! Print Imag(T)

!
! Calculate the energy derivative of the T_matrix
IF(EDeriv)THEN
   TE_Real=-SE_Real                                                                                     ! Calculate d(Real(T))/dE
   CALL Matrix_Out(TE_Real,NOpen,NOpen,'d(T_Real)/dE', 'd(Real Part of T-Matrix)/dE', Print_Flag)       ! Print d(Real(T))/dE
   TE_Imag=-SE_Imag                                                                                     ! Calculate d(Imag(T))/dE
   CALL Matrix_Out(TE_Imag,NOpen,NOpen,'d(T_Imag)/dE', 'd(Imaginary Part of T-Matrix)/dE', Print_Flag)  ! Print d(Imag(T))/dE
ENDIF

WRITE(Out_Unit,*)'Deallocate Temporary Storage in T_Matrix'
DEALLOCATE(Unit_Mat)

WRITE(Out_Unit,*)'Leaving:', ProcName

RETURN
ENDSUBROUTINE T_Matrix
