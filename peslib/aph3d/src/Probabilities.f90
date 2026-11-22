SUBROUTINE Probabilities (NOpen, S_Real, S_Imag, Prob, SE_Real, SE_Imag, ProbE, EDeriv)
USE Numeric_Kinds_Module
USE FileUnits_Asymptotic_Module
!
! Author: Gregory A. Parker, Department of Physics and Astronomy, University of Oklahoma
!
! Calculate the Probabilities and derivative of Probabilities with respect to Energy.
!
! Required Input <=====
!     NOpen           Size of ALL matrices (NOpen by NOpen)
!     EDeriv=True If ernergy derivatives of the K-Matrix exist
!     S_Real      Real part of the S-Matrix
!     S_Imag      Imaginary part of the S-Matrix
!     SE_Real     Energy derivative of the real part of the S-Matrix
!     SE_Imag     Energy derivative of the imaginary part of the S-Matrix
!
! On return  =====>
!     Prob        Probability matrix
!     ProbE       Energy derivative of the probability matrix
!
! This routine is called by:
!      Asymptotic
! This routine calls:
!      Matrix_Out
IMPLICIT NONE
CHARACTER(LEN=13), PARAMETER:: ProcName='Probabilities'
CHARACTER(LEN=6) Print_Flag
LOGICAL, INTENT(IN):: EDeriv
INTEGER, INTENT(IN):: NOpen
REAL(Kind=WP_Kind), INTENT(IN):: S_Real(NOpen,NOpen),  S_Imag(NOpen,NOpen)
REAL(Kind=WP_Kind), INTENT(IN):: SE_Real(NOpen,NOpen), SE_Imag(NOpen,NOpen)
REAL(Kind=WP_Kind), INTENT(OUT):: Prob(NOpen,NOpen), ProbE(NOpen,NOpen)

WRITE(Out_Unit,*)
WRITE(Out_Unit,*)
WRITE(Out_Unit,*)
WRITE(Out_Unit,*)'Entering:', ProcName
CALL PoptAsy(ProcName, Print_Flag)

!------------------------------------------------------------------------
! Calculate and print the probabilities.
!-----------------------------------------------------------------------
Prob=S_Real**2+S_Imag**2
CALL Matrix_Out(Prob,NOpen,NOpen,'Prob','Probability Matrix', Print_Flag)

!------------------------------------------------------------------------
! Calculate and print the energy derivative of the probabilities.
!-----------------------------------------------------------------------
IF(EDeriv)THEN
   ProbE=SE_Real*S_Real+S_Real*SE_Real+SE_Imag*S_Imag+S_Imag*SE_Imag
   CALL Matrix_Out(ProbE,NOpen,NOpen,'ProbE','Energy Derivative of Probability Matrix', Print_Flag)
ENDIF

WRITE(Out_Unit,*)'Leaving:', ProcName
RETURN
ENDSUBROUTINE Probabilities
