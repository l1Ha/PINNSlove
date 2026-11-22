SUBROUTINE Write_Vector2(avec, n, etot, ithmat, TypeOfMat, outunit)
USE Numeric_Kinds_Module
USE FileUnits_Asymptotic_Module
IMPLICIT NONE
INTEGER n
INTEGER ithmat, outunit
CHARACTER(LEN=20) TypeOfMat
REAL(Kind=WP_Kind) etot, avec(n)

WRITE(outunit)ithmat, TypeOfMat, n, etot
CALL writevec(avec, n, outunit)

RETURN
!
! This is reached at the end of the input file.
!
ENDSUBROUTINE Write_Vector2