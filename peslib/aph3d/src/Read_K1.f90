SUBROUTINE Read_K1(n, InUnit, OutUnit)
USE Numeric_Kinds_Module
USE ReadK_Module
IMPLICIT NONE
INTEGER n, OutUnit,InUnit
LOGICAL, PARAMETER:: DebugPrint=.False.
Character(LEN=20)Type1
CHARACTER(LEN=7), PARAMETER:: ProcName='Read_K1'
CHARACTER(LEN=6) Print_Flag

CALL PoptAsy(ProcName, Print_Flag)
READ(InUnit)ithmat, type1, n, etot
IF(DebugPrint)THEN
   WRITE(OutUnit,'(2x,A,I5,3A,I5,A,ES21.13)')"ithmat=",  ithmat, " type1=", TRIM(type1), " n=", n, " etot=", etot
   WRITE(*,'(2x,A,I5,3A,I5,A,ES21.13)')"ithmat=",  ithmat, " type1=", TRIM(type1), " n=", n, " etot=", etot
ENDIF


RETURN
ENDSUBROUTINE Read_K1