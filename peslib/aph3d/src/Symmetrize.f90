SUBROUTINE Symmetrize(NOpen, K_Matx, KE_Matx, EDeriv)
USE Numeric_Kinds_Module
USE FileUnits_Asymptotic_Module
USE Numbers_Module
IMPLICIT NONE
CHARACTER(LEN=10), PARAMETER:: ProcName='Symmetrize'
CHARACTER(LEN=6) Print_Flag

! Author: Gregory A. Parker, Department of Physics and Astronomy, University of Oklahoma
!
! Symmetrize the K-Matrix
!
! Required Input <=====
!     NOpen           Size of ALL matrices (NOpen by NOpen)
!     EDeriv=True If ernergy derivatives of the K-Matrix exist
!     K_Matx      K-Matrix
!     KE_Matx     Energy derivative of the K-Matrix
!
! On RETURN  =====>
!     K_Matx      K-Matrix (K)
!     KE_Matx     Energy derivative of the K-Matrix (K)
! This routine is called by:
!      Asymptotic
! This routine calls:
!      Matrix_Out
LOGICAL, INTENT(IN):: EDeriv
INTEGER, INTENT(IN):: NOpen
REAL(Kind=WP_Kind), INTENT(INOUT):: K_Matx(NOpen,NOpen), KE_Matx(NOpen,NOpen)

WRITE(Out_Unit,*)
WRITE(Out_Unit,*)
WRITE(Out_Unit,*)
WRITE(Out_Unit,*)'Entering:', ProcName
CALL PoptAsy(ProcName, Print_Flag)

!
!     Symmetrize the K-Matrix
K_Matx=Half*(K_Matx+Transpose(K_Matx))

!
!     Symmetrize the energy derivative of the K-Matrix
IF(EDeriv)THEN
   KE_Matx=Half*(KE_Matx+Transpose(KE_Matx))
ENDIF

!
! Print symmetrized K-Matrix
CALL Matrix_Out(K_Matx,NOpen,NOpen,  'K_Matrix',  'Symmetrized K_Matrix', Print_Flag)

!
! Print the symmetrized energy derivative of the K-Matrix
IF(EDeriv)THEN
   CALL Matrix_Out(KE_Matx,NOpen,NOpen,  'KE_Mat',  'Derivative of K_Matrix',  Print_Flag)
ENDIF

WRITE(Out_Unit,*)'Leaving:', ProcName
RETURN
ENDSUBROUTINE Symmetrize
