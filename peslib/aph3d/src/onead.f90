SUBROUTINE onead ( n, a)
USE FileUnits_Module
!-----------------------------------------------------------------------
!   This routine adds the identity matrix to matrix a.
!-----------------------------------------------------------------------
IMPLICIT NONE
!#include <Inc.unitno>
REAL(Kind=WP_Kind) a
INTEGER n, i
DIMENSION a(n,n)

!-----------------------------------------------------------------------
!  Add the Identity matrix.
!-----------------------------------------------------------------------
DO i = 1, n
    a(i,i) = a(i,i)+1.d0
ENDDO

RETURN
END
