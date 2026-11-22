SUBROUTINE VScale( n, a, ia, s )
USE Numeric_Kinds_Module
!
!-----------------------------------------------------------------------
!  this routine scales a(1), a(1+ia), ... a(1+n*ia) elements by s.
!-----------------------------------------------------------------------
IMPLICIT NONE
INTEGER n, ia ,i, ix
REAL(Kind=WP_Kind) a(n), s
!-----------------------------------------------------------------------
!  general case (ia not equal to 1).
!-----------------------------------------------------------------------
IF( ia /= 1 )THEN
   ix = 1
   DO i = 1, n
      a(ix) = a(ix)*s
      ix = ix + ia
   ENDDO
!-----------------------------------------------------------------------
!  special case for ia equal to 1.
!-----------------------------------------------------------------------
ELSE
   DO i = 1, n
      a(i) = a(i)*s
   ENDDO
ENDIF

RETURN
ENDSUBROUTINE VScale
