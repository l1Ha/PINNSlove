SUBROUTINE specadd(n, z, u)
  USE Numeric_Kinds_Module
  USE FileUnits_Module
  IMPLICIT NONE
!#include <Inc.unitno>
  INTEGER index, n, i, j
  REAL(Kind=WP_Kind) z(n,n), u(n*(n+1)/2)

  index = 0
  DO 20 i = 1, n
      DO 10 j = 1, i
         index = index + 1
         z(i,j) = z(i,j) + u(index)
         z(j,i) = z(i,j)
10      CONTINUE
20  CONTINUE
  RETURN
  END
