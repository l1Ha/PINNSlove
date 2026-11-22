SUBROUTINE MxOut_Unit(a, nrow, ncol, Print_Unit)
USE Numeric_Kinds_Module
USE FileUnits_Module
IMPLICIT NONE
CHARACTER(LEN=32), PARAMETER:: fmt1='(3x, "  row   ", 10(3x, i3, 6x))'
CHARACTER(LEN=35), PARAMETER:: fmt2='(1x, "col ", i3, 2x, 10000(es12.4))'
INTEGER nrow, ncol, irow, jcol, Print_Unit
CHARACTER(LEN=21), PARAMETER:: ProcName='MxOut_Unit'
CHARACTER(LEN=6) Print_Flag
REAL(Kind=WP_Kind) a(nrow,ncol)

CALL PoptAsy(ProcName, Print_Flag)
irow = 1
WRITE(Print_Unit,fmt1)(jcol,jcol=1,min(ncol,10))                     ! Write Column Labels
DO irow = 1, min(nrow,10)
   WRITE(Print_Unit,fmt2)irow, (a(irow,jcol),jcol=1,min(ncol,10))    ! Write Row label and Matrix Elements
ENDDO

RETURN
ENDSUBROUTINE MxOut_Unit
