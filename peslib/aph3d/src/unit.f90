SUBROUTINE Unit(n, UnitMat)
USE Numeric_Kinds_Module
IMPLICIT NONE
INTEGER n, i
REAL(Kind=WP_Kind) UnitMat(n,n)

UnitMat=0.d0
DO i=1,n
    UnitMat(i,i)=1.d0
ENDDO

RETURN
ENDSUBROUTINE Unit
