SUBROUTINE Write_Matrix2Old(amat, energies, chanl, elect, nvib, jrot,    &
                         lorb, n, etot, ithmat, TypeOfMat, jtot,      &
                         proj, approx, symm, outunit)
USE Numeric_Kinds_Module
USE FileUnits_Module
IMPLICIT NONE
INTEGER n, nstate, jtot, chanl(n), elect(n), nvib(n), jrot(n)
INTEGER lorb(n), ithmat, outunit
CHARACTER(LEN=20) TypeOfMat, approx, proj, symm
REAL(Kind=WP_Kind) etot, amat(n,n), energies(n)

WRITE(outunit)ithmat, TypeOfMat, n, etot, jtot, proj, approx, symm
CALL WriteMatOld (amat, energies, chanl, elect, nvib, jrot, lorb, n, outunit)

RETURN
!
! This is reached at the end of the input file.
!
ENDSUBROUTINE Write_Matrix2Old
