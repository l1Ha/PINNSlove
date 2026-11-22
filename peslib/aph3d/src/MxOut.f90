SUBROUTINE MxOut(a, nrow, ncol)
USE Numeric_Kinds_Module
USE FileUnits_Module
IMPLICIT NONE
CHARACTER(LEN=23),     PARAMETER:: fmt0='(A,I4,A,I4,A,I5,A,I5,A)'
CHARACTER(LEN=32),     PARAMETER:: fmt1='(3x, "  row   ", 10(3x, i3, 6x))'
CHARACTER(LEN=32),     PARAMETER:: fmt2='(1x, "col ", i3, 2x, 10(es12.4))'
INTEGER(Kind=IW_Kind), PARAMETER:: MatPrint=10
CHARACTER(LEN=21), PARAMETER:: ProcName='MxOut'
INTEGER nrow, ncol, irow, jcol
CHARACTER(LEN=6) Print_Flag
REAL(Kind=WP_Kind) a(nrow,ncol)

CALL PoptAsy(ProcName, Print_Flag)
irow = 1
WRITE(Out_Unit,fmt0) "Printing a ", MatPrint, " by", MatPrint, "Block of a", nrow, " by", ncol, " Matrix"
WRITE(Out_Unit,fmt1)(jcol,jcol=1,min(ncol,MatPrint))                     ! Write Column Labels
DO irow = 1, min(nrow,MatPrint)
   WRITE(Out_Unit,fmt2)irow, (a(irow,jcol),jcol=1,min(ncol,MatPrint))    ! Write Row label and Matrix Elements
ENDDO

RETURN
ENDSUBROUTINE MxOut
