SUBROUTINE vsets ( n, a, ia, s )
USE Numeric_Kinds_Module
USE FileUnits_Module
!

!  p u r p o s e   o f    s u b r o u t i n e
!-----------------------------------------------------------------------
! this routine sets a(1), a(1+ia), ... a(1+n*ia) elements equal to s.
!-----------------------------------------------------------------------


!  i n p u t    a r g u m e n t s
!  n
!  a
!  ia
!  s

IMPLICIT NONE
INTEGER n, ia ,i, ix

!  r e a l s
REAL(Kind=WP_Kind) a, s

!  d i m e n s i o n s

DIMENSION a(n)

!-----------------------------------------------------------------------
!  general case (ia not equal to 1).
!-----------------------------------------------------------------------
IF( ia /= 1 )THEN
   ix = 1
   DO 1 i = 1, n
      a(ix) = s
      ix = ix + ia
 1       CONTINUE
!-----------------------------------------------------------------------
!  special case for ia equal to 1.
!-----------------------------------------------------------------------
ELSE
   DO 2 i = 1, n
      a(i) = s
 2       CONTINUE
ENDIF

RETURN
END
