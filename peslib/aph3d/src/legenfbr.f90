SUBROUTINE legenfbr (norder, pn, x)

!            P U R P O S E   O F    S U B R O U T I N E
!   This routine calculates normalized Legendre polynomials at x.


!            I N P U T    A R G U M E N T S
!   norder n order of the largest Legendre polynomial desired.
!   x      argument of the polynomial in the range -1<=x<=1.

!            O U T P U T    A R G U M E N T S
!   pn     array containing the Legendre polynomials evaluated at x.
!   dpndx  array containing first derivatives of the pn.


USE FileUnits_Module
USE Numeric_Kinds_Module
USE Numbers_Module
IMPLICIT NONE
INTEGER  n, norder
REAL(Kind=WP_Kind)pn(0:norder), x, sqr      !, dpndx(0:norder)
!-----------------------------------------------------------------------
! Determine printing options.
!-----------------------------------------------------------------------
LOGICAL little, medium, full
INTEGER ithcall, ithsub
DATA ithcall /0/, ithsub /0/
DATA little /.false./, medium /.false./, full /.false./
CALL popt ('legendre', little, medium, full, ithcall, ithsub)
!-----------------------------------------------------------------------
! Check the input arguments for consistency.
!-----------------------------------------------------------------------
IF(norder<0.or.ABS(x)>one)THEN
   WRITE(Out_Unit,*)'Incorrect arguments to legendre'
   WRITE(Out_Unit,*)'norder = ',norder,'  x =',x
   WRITE(*,*)'Incorrect arguments to legendre'
   WRITE(*,*)'norder = ',norder,'  x =',x
   STOP 'Legenfbr'
ENDIF

!-----------------------------------------------------------------------
!  Evaluate the Legendre polynomials. This is the standard recurrsion
!  relation rewritten for numerical stability.
!-----------------------------------------------------------------------
pn(0) = one
!     dpndx(0)=zero
IF(norder>0)THEN
   pn(1) = x
!        dpndx(1)=one
ENDIF
IF(norder>1)THEN
   DO n = 2, norder
    pn(n) = x*pn(n-1)-pn(n-2)+x*pn(n-1)-(x*pn(n-1)-pn(n-2))/n
!*         dpndx(n)=n*(x*pn(n)-pn(n-1))/(x**2-one)
   ENDDO
ENDIF

!  -------------------------------------------------------------------
!  normalization
!  -------------------------------------------------------------------
DO n=0,norder
   sqr=sqrt(n+half)
   pn(n)=sqr*pn(n)
!*     dpndx(n)=sqr*dpndx(n)
ENDDO

!-----------------------------------------------------------------------
! IF desired WRITE out the Legendre polynomials.
!-----------------------------------------------------------------------
IF(medium)THEN
   WRITE(Out_Unit,*)'routine Legendre'
   WRITE(Out_Unit,*)'maximum order = ',norder
ENDIF
IF(full)THEN
   WRITE(Out_Unit,*)'Legendre polynomials evaluated at x = ',x
   DO n = 0, norder
      WRITE(Out_Unit,*)'n= ',n,'  Pn(x) = ',pn(n)
   ENDDO
ENDIF

RETURN
END
