SUBROUTINE fill_mat(a, n)
USE Numeric_Kinds_Module
USE FileUnits_Module
IMPLICIT NONE
!#include <Inc.unitno>
INTEGER i, j, n
REAL(Kind=WP_Kind) a(n,n)

DO i = 1, n
   DO j = 1, i
      a(j,i) = a(i,j)
   ENDDO
ENDDO

RETURN
END
