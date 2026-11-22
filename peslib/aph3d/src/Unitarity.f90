SUBROUTINE Unitarity(NOpen, S_Real, S_Imag, SE_Real, SE_Imag, EDeriv)
USE Numeric_Kinds_Module
USE FileUnits_Asymptotic_Module
USE matrixlabels_Module      ! Use matrixlabels module for printing quantum labels
IMPLICIT NONE
CHARACTER(LEN=9), PARAMETER:: ProcName='Unitarity'      ! Procedure name
CHARACTER(LEN=6) Print_Flag                             ! Print Flag
LOGICAL DBug/.False./                        ! Debugging parameter

! Author: Gregory A. Parker, Department of Physics and Astronomy, University of Oklahoma
!
! Determines unitarity of the S-Matrix
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
!     Nothing new is returned
!
! This routine is called by:
!      Asymptotic
! This routine calls:
!      LgDiff
!
LOGICAL, INTENT(IN):: EDeriv                                             ! True if Energy derivatives were calculated
INTEGER, INTENT(IN):: NOpen                                              ! Number of States
INTEGER IState, JState                                                   ! Local indexing variables
REAL(Kind=WP_Kind), INTENT(IN)::  S_Real(NOpen,NOpen),  S_Imag(NOpen,NOpen)  ! S-Matrix
REAL(Kind=WP_Kind), INTENT(IN):: SE_Real(NOpen,NOpen), SE_Imag(NOpen,NOpen)  ! Energy derivative of S-Matrix
REAL(Kind=WP_Kind) BigDiff
REAL(Kind=WP_Kind), ALLOCATABLE:: UniCheck(:), UniCheckE(:)

WRITE(Out_Unit,*)
WRITE(Out_Unit,*)
WRITE(Out_Unit,*)
WRITE(Out_Unit,*)'Entering:', ProcName
CALL PoptAsy(ProcName, Print_Flag)
IF(Print_Flag=='Full')DBug=.True.

ALLOCATE(UniCheck(NOpen),UniCheckE(NOpen))
LabelRows=.False.

!
! Determine row unitarities
!
DO IState = 1, NOpen
   UniCheck(IState) = 0.d0
   IF(EDeriv)UniCheckE(IState)=0.d0
   DO JState = 1, NOpen
      UniCheck(IState) = UniCheck(IState) + S_Real(IState,JState)**2 &
                                          + S_Imag(IState,JState)**2
      IF(EDeriv)THEN
         UniCheckE(IState)=UniCheckE(IState)+S_Real(IState,JState)*SE_Real(IState,JState)  &
                                            +S_Imag(IState,JState)*SE_Imag(IState,JState)
      ENDIF
   ENDDO
ENDDO
IF(EDeriv)UniCheckE=2.d0*UniCheckE

!
! Determine largest percent non-unitarity for each row and print results
IF(Print_Flag/='None')THEN
   WRITE(Out_Unit,*)'Row Checks'
   CALL LgDiff(NOpen, UniCheck, BigDiff)
   IF(DBug.or.BigDiff>1.d-6)THEN
      CALL Matrix_Out(UniCheck, 1, NOpen, 'UniCheck', 'Row Unitarities', Print_Flag)
      UniCheck=UniCheck-1.d0
      CALL Matrix_Out(UniCheck, 1, NOpen, 'UniCheck', 'Differences from Unity of Rows', Print_Flag)
   ENDIF
   WRITE(Out_Unit,*)'Largest percent nonunitarity=', BigDiff


!
! Energy derivative of unitarities should be zero
! Calculate largest percent row error
   IF(EDeriv)THEN
     CALL Matrix_Out(UniCheckE, 1, NOpen, 'UniCheckE', 'Row D[sum,E] check (Should be Zero)', Print_Flag)
     BigDiff=100.d0*MaxVal(ABS(UniCheckE))
     WRITE(Out_Unit,*)'Largest error=', BigDiff
   ENDIF
ENDIF

LabelRows=.True.
LabelColumns=.False.
!
! Determine Column unitarities
DO JState = 1, NOpen
  UniCheck(JState) = 0.d0
  IF(EDeriv)UniCheckE(JState)=0.d0
  DO IState = 1, NOpen
     UniCheck(JState) = UniCheck(JState) + S_Real(IState,JState)**2 + S_Imag(IState,JState)**2
     IF(EDeriv)THEN
       UniCheckE(JState)=UniCheckE(JState)+S_Real(IState,JState)*SE_Real(IState,JState)+S_Imag(IState,JState)*SE_Imag(IState,JState)
     ENDIF
  ENDDO
ENDDO
IF(EDeriv)UniCheckE=2.d0*UniCheckE

!
! Determine largest percent non-unitarity for each column and print results
IF(Print_Flag/='None')THEN
   WRITE(Out_Unit,*)'Column Checks'
   CALL LgDiff(NOpen, UniCheck, BigDiff)
   IF(DBug.or.BigDiff>1.d-6)THEN
      CALL Matrix_Out(UniCheck, NOpen, 1, 'UniCheck', 'Column Unitarities', Print_Flag)
      UniCheck=UniCheck-1.d0
      CALL Matrix_Out(UniCheck, NOpen, 1, 'UniCheck', 'Differences from Unity of Columns', Print_Flag)
      WRITE(Out_Unit,*)'Largest percent nonunitarity=', BigDiff
   ENDIF

!
! Energy derivative of unitarities should be zero
! Calculate largest percent column error
   IF(EDeriv)THEN
      CALL Matrix_Out(UniCheckE, NOpen, 1, 'UniCheckE', 'Column D[sum,E] check (Should be Zero)', Print_Flag)
      BigDiff=100.0d0*MaxVal(ABS(UniCheckE))
      WRITE(Out_Unit,*)'Largest error=', BigDiff
   ENDIF
ENDIF

LabelRows=.True.
LabelColumns=.True.
WRITE(Out_Unit,*)'Deallocate Temporary Storage in Unitarity'
DEALLOCATE(UniCheck,UniCheckE)

WRITE(Out_Unit,*)'Leaving:', ProcName
RETURN
ENDSUBROUTINE Unitarity
