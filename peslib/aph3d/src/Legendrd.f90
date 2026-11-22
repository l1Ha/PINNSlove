SUBROUTINE Legendrd (norder, m, pnm, dpnmdx, x, nfac)
!

!  p u r p o s e   o f    s u b r o u t i n e
!  calculates, for a given x and m, all the normalized associated
!  legendre polynomials pnm from n=m up to n=norder.
!  also calculates their first derivatives.

!  i n p u t    a r g u m e n t s
!  norder n order of the largest legendre polynomial desired.
!  x      argument of the polynomial in the range -1<=x<=1.

!  o u t p u t    a r g u m e n t s
!  pnm     array of normed, assoc legendre polynomials evaluated at x.
!  dpnmdx  array containing first derivatives of the pnm.

!  <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

USE Numeric_Kinds_Module
USE Parms_Module
USE FileUnits_Module
USE Numbers_Module
IMPLICIT NONE
SAVE
INTEGER  n, m, norder, nthcall, i
REAL(Kind=WP_Kind) x, sqr
REAL(Kind=WP_Kind) fact, omx2, somx2, abs
REAL(Kind=WP_Kind) pnm(0:mxl), dpnmdx(0:mxl), nfac(0:2*mxl)
INTRINSIC abs, sqrt
REAL(Kind=WP_Kind), PARAMETER:: eps=1d-8

!-----------------------------------------------------------------------
!  determine printing options.
!-----------------------------------------------------------------------
LOGICAL little, medium, full
INTEGER ithcall, ithsub
DATA nthcall /0/
DATA ithcall /0/, ithsub /0/
DATA little /.false./, medium /.false./, full /.false./
CALL popt ('legendre', little, medium, full, ithcall, ithsub)
!-----------------------------------------------------------------------
!  check the input arguments for consistency.
!-----------------------------------------------------------------------
IF(m < 0 .or. m > norder .or. ABS(x) > one+eps)THEN
   WRITE(Out_Unit,*)'incorrect arguments to legendre'
   WRITE(Out_Unit,*)'norder = ', norder, ' m = ', m, ' x =', x
   STOP 'legendre'
ENDIF
IF(ABS(x) == one)THEN
!  WRITE(Out_Unit,*)' ***warning***  legendre called with ABS(x)=1'
!  WRITE(Out_Unit,*)' pnm will be correct, but dpnmdx may be incorrect.'
ENDIF

!IF(nthcall > 0) GOTO 4
!  ---------------------------------------------------------------
!  calculate factorials for norm on first pass
!  ---------------------------------------------------------------
nfac(0)=one
DO n=1, 2*mxl
   nfac(n)=n*nfac(n-1)
ENDDO
nthcall=1

4 CONTINUE
!-----------------------------------------------------------------------
!  evaluate the legendre polynomials and their first derivative.
!  pnm recursion formula is from numerical recipes (1986), p. 181.
!  the formula for dpnmdx is gradshteyn-ryzhik (1965), 8.733.1b.
!  note.  dpnmdx is not exact at x=1 or -1, but that is not used by
!  gaussian quadratures.
!  first compute pnm(n=m) and derivative.  DO m=0 first, THEN m>0.
!-----------------------------------------------------------------------
pnm(m) = one
dpnmdx(m)=zero
omx2=(one-x)*(one+x)
somx2=sqrt(omx2)
IF(omx2 == zero)omx2=1.d-100
IF(m > 0)THEN
   fact=one
   DO i=1, m
      pnm(m)=-pnm(m)*fact*somx2
      fact=fact+two
   ENDDO
   dpnmdx(m)=-m*x*pnm(m)/omx2
ENDIF
IF(norder > m)THEN
!  ---------------------------------------------------------------
!  compute pnm(m+1) and derivative.
!  ---------------------------------------------------------------
   pnm(m+1)=x*(2*m+1)*pnm(m)
   dpnmdx(m+1)=(-(m+1)*x*pnm(m+1)+(2*m+1)*pnm(m))/omx2
ENDIF
IF(norder > m+1)THEN
!  ---------------------------------------------------------------
!  compute pnm(n) and derivative for n>m+1
!  ---------------------------------------------------------------
   DO n = m+2, norder
      pnm(n)=(x*(2*n-1)*pnm(n-1)-(n+m-1)*pnm(n-2))/(n-m)
      dpnmdx(n)=(-n*x*pnm(n)+(n+m)*pnm(n-1))/omx2
   ENDDO
ENDIF

!  -------------------------------------------------------------------
!  at this point, pnm is the standard unnormalized assoc. legendre poly.
!  ---------------------------------------------------------------
!  normalization
!  -------------------------------------------------------------------
DO n=m, norder
   sqr=sqrt((n+half)*nfac(n-m)/nfac(n+m))
   pnm(n)=sqr*pnm(n)
   dpnmdx(n)=sqr*dpnmdx(n)
ENDDO

!-----------------------------------------------------------------------
!  IF desired WRITE out the legendre polynomials.
!-----------------------------------------------------------------------
IF(medium)THEN
   WRITE(Out_Unit,*)'routine legendre'
   WRITE(Out_Unit,*)'maximum order = ', norder
ENDIF
IF(full)THEN
   WRITE(Out_Unit,*)'legendre polynomials evaluated at x = ', x
   DO n = m, norder
      WRITE(Out_Unit,*)'n= ', n, ' m=', m, '  pnm(x) = ', pnm(n)
   ENDDO
ENDIF

RETURN
ENDSUBROUTINE Legendrd
