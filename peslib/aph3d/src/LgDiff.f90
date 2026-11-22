SUBROUTINE LgDiff(NOpen, UniCheck, BigDiff)
USE Numeric_Kinds_Module
USE FileUnits_Asymptotic_Module
IMPLICIT NONE
CHARACTER(LEN=21), PARAMETER:: ProcName='LgDiff'  ! Procedure name
CHARACTER(LEN=6) Print_Flag
LOGICAL, PARAMETER:: DBug=.false.                ! Debugging parameter

!
! Author: Gregory A. Parker, Department of Physics and Astronomy, University of Oklahoma
! Calculates largest percent non-unitarity.
!
!     Required Input <=====
!     NOpen           Size of ALL matrices (NOpen by NOpen)
!     UniCheck        An array of unitarities
!
!     On RETURN  =====>
!     BigDiff         Largest percent non-unitarity
!
! This routine is called by:
!      Unitarity
INTEGER, INTENT(IN):: NOpen                    ! Number of States
INTEGER IState                                 ! Local index
REAL(Kind=WP_Kind), INTENT(IN):: UniCheck(NOpen)   ! Unitarity
REAL(Kind=WP_Kind), INTENT(OUT)::BigDiff           ! Largest percent non-unitarity
REAL(Kind=WP_Kind) Diff                           ! Local variable

CALL PoptAsy(ProcName, Print_Flag)
IF(DBug)THEN
   WRITE(Out_Unit,*)
   WRITE(Out_Unit,*)
   WRITE(Out_Unit,*)
   WRITE(Out_Unit,*)'Entering:', ProcName
ENDIF

!
! Calculate largest percent non-unitarity
BigDiff = 0.d0
DO  IState = 1, NOpen
   Diff = 100.0d0*ABS(UniCheck(IState)-1.d0)
   IF(BigDiff<Diff) BigDiff = Diff
ENDDO

IF(DBug)THEN
   WRITE(Out_Unit,*)'Leaving:', ProcName
ENDIF

RETURN
ENDSUBROUTINE LgDiff
