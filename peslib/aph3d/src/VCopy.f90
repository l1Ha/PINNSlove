SUBROUTINE VCopy( n, a, b)
USE Numeric_Kinds_Module
!
!-----------------------------------------------------------------------
!  this routine copies the vector b into the vector a
!-----------------------------------------------------------------------
IMPLICIT NONE
REAL(Kind=WP_Kind) a, b
INTEGER n, i
DIMENSION a(n), b(n)

!-----------------------------------------------------------------------
!  copy the vector.
!-----------------------------------------------------------------------
DO i = 1, n
   a(i) = b(i)
ENDDO

RETURN
ENDSUBROUTINE VCopy
