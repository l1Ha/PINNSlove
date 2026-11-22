SUBROUTINE dhep (h, x, n)
USE Numeric_Kinds_Module
USE FileUnits_Module
INTEGER I, N
REAL(Kind=WP_Kind) h(n+1), x
!-------------------------------------------------------------------
!  written by g. a. parker
!  recuresion to calculate the hermite polynomials h(i+1) for
!  i=0 to i=n at x.
!-------------------------------------------------------------------
REAL(Kind=WP_Kind) hnp1, hn, hnm1, sqr2, sqr2dx, dx
REAL(Kind=WP_Kind) pi4i
DATA sqr2 /1.4142135623730950488016887d0/
DATA pi4i /7.511255444649424828587030d-1/
dx=x
sqr2dx=sqr2*dx
hnm1=pi4i
hn=sqr2dx*pi4i
h(1)=hnm1
h(2)=hn
IF(n >= 2) GOTO 10
RETURN
!-------------------------------------------------------------------
!  start recursion
!-------------------------------------------------------------------
 10   CONTINUE
 DO i=2,n
   hnp1=(sqr2dx*hn-sqrt( REAL(i-1,WP_Kind) )*hnm1)/sqrt( REAL(i,WP_Kind) )
   hnm1=hn
   h(i+1)=hnp1
   hn=hnp1
ENDDO
RETURN
!----------------***END-dhep***---------------------------------------
END
