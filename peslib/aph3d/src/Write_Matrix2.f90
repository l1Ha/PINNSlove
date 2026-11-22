SUBROUTINE Write_Matrix2(amat, n, etot, ithmat, TypeOfMat, outunit)
USE FileUnits_Asymptotic_Module
IMPLICIT NONE
INTEGER n
INTEGER ithmat, outunit
CHARACTER(LEN=20) TypeOfMat
CHARACTER(LEN=13), PARAMETER:: ProcName='Write_Matrix2'
CHARACTER(LEN=6) Print_Flag
REAL(Kind=WP_Kind) etot, amat(n,n)

CALL PoptAsy(ProcName, Print_Flag)
WRITE(outunit)ithmat, TypeOfMat, n, etot
CALL WriteMat(amat, n, outunit)

RETURN
!
! This is reached at the end of the input file.
!
ENDSUBROUTINE Write_Matrix2
