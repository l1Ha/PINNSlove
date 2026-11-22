SUBROUTINE Legendre (norder, m, pnm, x)
USE Numeric_Kinds_Module

!  P U R P O S E   O F    S U B R O U T I N E
!  Calculates, for a given x and m, all the normalized associated
!  Legendre polynomials Pnm from n=m up to n=norder.
!  R T Pack, Dec. 1988.  Another version written at the same time for
!  aph analytic basis calcs, also calcs the first derivative.  This
!  version does not.


!  I N P U T    A R G U M E N T S
!  norder n order of the largest Legendre polynomial desired.
!  x      argument of the polynomial in the range -1<=x<=1.

!  O U T P U T    A R G U M E N T S
!  pnm     array of normed, assoc Legendre polynomials evaluated at x.
USE FileUnits_Module
USE Narran_Module
USE Parms_Module

IMPLICIT NONE
SAVE
INTEGER  n, m, norder, nthcall, i
REAL(Kind=WP_Kind) x, sqr, zero, half, one, two, fact, omx2, somx2, eps
REAL(Kind=WP_Kind) pnm(0:mxl)
REAL(Kind=WP_Kind), ALLOCATABLE:: nfac(:)
INTRINSIC abs, sqrt
EXTERNAL popt
PARAMETER (eps=1d-8, zero=0.0d0, half=0.5d0, one=1.0d0,two=2.0d0)

DATA nthcall /0/

!-----------------------------------------------------------------------
!  Determine printing options.
!-----------------------------------------------------------------------
LOGICAL little, medium, full
INTEGER ithcall, ithsub
DATA ithcall /0/, ithsub /0/
DATA little /.false./, medium /.false./, full /.false./
CALL popt ('legendre', little, medium, full, ithcall, ithsub)
!-----------------------------------------------------------------------
!  Check the input arguments for consistency.
!-----------------------------------------------------------------------
IF(m < 0 .or. m > norder .or. ABS(x) > one+eps)THEN
   WRITE(Out_Unit,*)'Incorrect arguments to legendre'
   WRITE(Out_Unit,*)'norder = ', norder, ' m = ', m, ' x =', x
   WRITE(*,*)'Incorrect arguments to legendre'
   WRITE(*,*)'norder = ', norder, ' m = ', m, ' x =', x
   STOP 'Legendre'
ENDIF
IF(.NOT.allocated(nfac))ALLOCATE(nfac(0:2*mxl))
!  ---------------------------------------------------------------
!  calculate factorials for norm on first pass
!  ---------------------------------------------------------------
IF(nthcall == 0)THEN
   IF(.NOT.ALLOCATED(nfac))ALLOCATE(NFac(0:200))
   nfac(0)=one
   DO n=1, 2*mxl
      nfac(n)=n*nfac(n-1)
   ENDDO
   nthcall=1
ENDIF
!-----------------------------------------------------------------------
!  evaluate the legendre polynomials.
!  pnm recursion formula is from numerical recipes (1986), p. 181.
!  DO m=0 first, THEN m>0.
!  get pmm first THEN recur upward from there in n.
!-----------------------------------------------------------------------
pnm(m) = one
omx2=(one-x)*(one+x)
somx2=sqrt(omx2)
IF(m > 0)THEN
   fact=one
   DO 11 i=1, m
      pnm(m)=-pnm(m)*fact*somx2
      fact=fact+two
 11      CONTINUE
ENDIF
IF(norder > m)THEN
!  ---------------------------------------------------------------
!  compute pnm(m+1) .
!  ---------------------------------------------------------------
   pnm(m+1)=x*(2*m+1)*pnm(m)
ENDIF
IF(norder > m+1)THEN
!  ---------------------------------------------------------------
!  compute pnm(m+1) .
!  ---------------------------------------------------------------
   pnm(m+1)=x*(2*m+1)*pnm(m)
ENDIF
IF(norder > m+1)THEN
!  ---------------------------------------------------------------
!  compute pnm(n)  for n>m+1
!  ---------------------------------------------------------------
   DO 1 n = m+2, norder
      pnm(n)=(x*(2*n-1)*pnm(n-1)-(n+m-1)*pnm(n-2))/(n-m)
    1    CONTINUE
ENDIF
!  -------------------------------------------------------------------
!  at this point, pnm is the standard unnormalized assoc. legendre poly.
!  ---------------------------------------------------------------
!  normalization
!  -------------------------------------------------------------------
DO 2 n=m, norder
   sqr=sqrt((n+half)*nfac(n-m)/nfac(n+m))
   pnm(n)=sqr*pnm(n)
 2    CONTINUE

!-----------------------------------------------------------------------
!  IF desired WRITE out the Legendre polynomials.
!-----------------------------------------------------------------------
IF(medium)THEN
   WRITE(Out_Unit,*)'routine Legendre'
   WRITE(Out_Unit,*)'maximum order = ', norder
ENDIF
IF(full)THEN
   WRITE(Out_Unit,*)'Legendre polynomials evaluated at x = ', x
   DO 5 n = m, norder
      WRITE(Out_Unit,*)'n= ', n, ' m=', m, '  Pnm(x) = ', pnm(n)
    5    CONTINUE
ENDIF

RETURN
ENDSUBROUTINE Legendre
