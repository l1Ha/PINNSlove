SUBROUTINE Unit_Matrix(NOpen, Unit_Mat)
USE Numeric_Kinds_Module
USE FileUnits_Asymptotic_Module
IMPLICIT NONE
CHARACTER(LEN=21), PARAMETER:: ProcName='Unit_Matrix'    ! Procedure name
LOGICAL, PARAMETER:: DBug=.false.                        ! Debugging parameter
CHARACTER(LEN=6) Print_Flag

! Author: Gregory A. Parker, Department of Physics and Astronomy, University of Oklahoma
!
! Creates an NOpen-by-NOpen Unit Matrix.
!
!     Required Input <=====
!     NOpen           Size of ALL matrices (NOpen by NOpen)
!
!     On return  =====>
!     Unit_Mat     An NOpen-by-NOpen unit matrix
!
! This routine is called by:
!      S_Matrix
!      T_Matrix
! This routine calls:
!      Matrix_Out
INTEGER, INTENT(IN):: NOpen                             ! Number of States
INTEGER IState                                           ! Local Index
REAL(Kind=WP_Kind), INTENT(OUT):: Unit_Mat(NOpen,NOpen)    ! Unit Matrix

IF(DBug)THEN
   WRITE(Out_Unit,*)
   WRITE(Out_Unit,*)
   WRITE(Out_Unit,*)
   WRITE(Out_Unit,*)'Entering:', ProcName
ENDIF

CALL PoptAsy(ProcName, Print_Flag)

!
! Create Unit Matrix
Unit_Mat=0.d0
DO IState=1,NOpen
   Unit_Mat(IState,IState)=1.d0
ENDDO

IF(DBug)THEN
   CALL Matrix_Out(Unit_Mat, NOpen, NOpen, 'Unit_Mat', 'Unit Matrix',  Print_Flag)
   WRITE(Out_Unit,*)'Leaving:', ProcName
ENDIF

RETURN
ENDSUBROUTINE Unit_Matrix
