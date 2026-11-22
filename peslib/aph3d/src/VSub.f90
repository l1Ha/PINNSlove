SUBROUTINE VSub( n, a, b)
USE Numeric_Kinds_Module
!
!-----------------------------------------------------------------------
!  this routine substracts vector b from vector a
!  and stores the result in a.
!-----------------------------------------------------------------------
IMPLICIT NONE
SAVE
REAL(Kind=WP_Kind) a, b
INTEGER n, i
DIMENSION a(n), b(n)

!-----------------------------------------------------------------------
!  substract the vectors.
!-----------------------------------------------------------------------
DO i = 1, n
   a(i) = a(i)-b(i)
ENDDO

RETURN
ENDSUBROUTINE VSub
