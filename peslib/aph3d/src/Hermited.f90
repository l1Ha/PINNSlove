SUBROUTINE Hermited (norder, hn, dhndx, x)
!

!  p u r p o s e   o f    s u b r o u t i n e
!  this routine calculates the hermite polynomials at x.


!  i n p u t    a r g u m e n t s
!  norder   maximum order of the hermite polynomial desired.
!  x        argument of the hermite polynomial.

!  o u t p u t    a r g u m e n t s
!  hn       array containing the hermite polynomials evaluated at x.
!  dhndx      array containing fcn needed for derivative of sho fcn.

!  <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>


USE Numeric_Kinds_Module
USE FileUnits_Module
USE Numbers_Module
IMPLICIT NONE
INTEGER  n, norder
REAL(Kind=WP_Kind) x
REAL(Kind=WP_Kind)  hn(0:norder), dhndx(0:norder)
EXTERNAL popt
!-----------------------------------------------------------------------
!  determine printing options.
!-----------------------------------------------------------------------
LOGICAL little, medium, full
INTEGER ithcall, ithsub
DATA ithcall /0/, ithsub /0/
DATA little /.false./, medium /.false./, full /.false./
CALL popt ('hermited', little, medium, full, ithcall, ithsub)
!-----------------------------------------------------------------------
!  check for input argument consistency.
!-----------------------------------------------------------------------
IF(norder < 0)THEN
   WRITE(Out_Unit,*)'incorrect arguments to hermited'
   WRITE(Out_Unit,*)'norder = ', norder, '  x =', x
   STOP 'hermited'
ENDIF
!-----------------------------------------------------------------------
!  evaluate the hermite polynomials using the standard recurrsion
!  relation.
!-----------------------------------------------------------------------
hn(0) = one
dhndx(0)=-x
IF(norder > 0)THEN
   hn(1) = two*x
   dhndx(1)=two-x*hn(1)
ENDIF
IF(norder > 1)THEN
   DO n = 2, norder
      hn(n) = two*(x*hn(n-1)-(n-1)*hn(n-2))
      dhndx(n)=two*n*hn(n-1) - x*hn(n)
   ENDDO
ENDIF

!-----------------------------------------------------------------------
!  IF desired WRITE out the hermite polynomials.
!-----------------------------------------------------------------------
IF(medium)THEN
   WRITE(Out_Unit,*)'hermite polynomials evaluated at x = ', x
   WRITE(Out_Unit,*)'maximum order = ', norder
ENDIF
IF(full)THEN
   DO n = 0, norder
      WRITE(Out_Unit,*)'n= ', n, '  hn(x) = ', hn(n)
   ENDDO
ENDIF
!
RETURN
ENDSUBROUTINE Hermited
