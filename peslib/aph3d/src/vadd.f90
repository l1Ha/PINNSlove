SUBROUTINE vadd ( n, a, b)
USE FileUnits_Module
!
!-----------------------------------------------------------------------
!  this routine adds vector a to vector b and stores the result in a.
!-----------------------------------------------------------------------
IMPLICIT NONE
REAL(Kind=WP_Kind) a, b
INTEGER n, i
DIMENSION a(n), b(n)

!-----------------------------------------------------------------------
!  add the vectors.
!-----------------------------------------------------------------------
DO i = 1, n
   a(i) = a(i)+b(i)
ENDDO

RETURN
ENDSUBROUTINE VAdd
