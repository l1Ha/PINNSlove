FUNCTION IndexR (String)
IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
! >>> IndexR -- returns the index of the last non-blank CHARACTER, or a
! >>>             if the entire string is blank
!
! *** Calling Sequence:         IPos    = IndexR ( String )
!
!-----------------------------------------------------------------------
!
!
INTEGER :: IndexR
CHARACTER(LEN=*) :: String
!
! * * * * *                             !       loop index
INTEGER :: i
!
!----------------------- execution begins here -------------------------
!
DO i = len (String), 1, - 1
   IF(String (i:i) /=' ')THEN
      IndexR = i
      RETURN
   ENDIF
ENDDO
!
IndexR = 0
!
RETURN
ENDFUNCTION IndexR
