SUBROUTINE vseti ( n, a, ia, s )
USE Numeric_Kinds_Module
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

!  o u t p u t    a r g u m e n t s

! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

USE FileUnits_Module
IMPLICIT NONE
INTEGER n, ia ,i, ix
INTEGER a(n), s

!-----------------------------------------------------------------------
!  general case (ia not equal to 1).
!-----------------------------------------------------------------------
IF( ia .ne. 1 )THEN
   ix = 1
   DO i = 1, n
      a(ix) = s
      ix = ix + ia
   ENDDO
!-----------------------------------------------------------------------
!  special case for ia equal to 1.
!-----------------------------------------------------------------------
ELSE
   DO i = 1, n
      a(i) = s
   ENDDO
ENDIF

RETURN
END
