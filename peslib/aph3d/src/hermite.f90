SUBROUTINE hermite (norder, hn, x)
USE Numeric_Kinds_Module
!

!  P U R P O S E   O F    S U B R O U T I N E
!  This routine calculates the Hermite polynomials at x.


!  I N P U T    A R G U M E N T S
!  norder   maximum order of the Hermite polynomial desired.
!  x        argument of the Hermite polynomial.

!  O U T P U T    A R G U M E N T S
!  hn       array containing the Hermite polynomials evaluated at x.
!  dhndx      array containing fcn needed for derivative of sho fcn.

 USE FileUnits_Module

IMPLICIT NONE


!  I N T E G E R S
INTEGER  n, norder

!  R E A L S
REAL(Kind=WP_Kind) x, hn, dhndx, one, two
PARAMETER (one=1.0d0, two=2.0d0)

!  D I M E N S I O N S
DIMENSION hn(0:norder)


!  E X T E R N A L S
EXTERNAL popt
!-----------------------------------------------------------------------
!  Determine printing options.
!-----------------------------------------------------------------------
LOGICAL little, medium, full
INTEGER ithcall, ithsub
DATA ithcall /0/, ithsub /0/
DATA little /.false./, medium /.false./, full /.false./
CALL popt ('hermite ', little, medium, full, ithcall, ithsub)
!-----------------------------------------------------------------------
!  Check for input argument consistency.
!-----------------------------------------------------------------------
IF(norder < 0)THEN
   WRITE(Out_Unit,*)'Incorrect arguments to hermite'
   WRITE(Out_Unit,*)'norder = ', norder, '  x =', x
   STOP 'hermite'
ENDIF
!-----------------------------------------------------------------------
!  Evaluate the Hermite polynomials using the standard recurrsion
!  relation.
!-----------------------------------------------------------------------
hn(0) = one
IF(norder > 0)THEN
   hn(1) = two*x
ENDIF
IF(norder > 1)THEN
   DO 1 n = 2, norder
      hn(n) = two*(x*hn(n-1)-(n-1)*hn(n-2))
    1    CONTINUE
ENDIF

!-----------------------------------------------------------------------
!  IF desired WRITE out the Hermite polynomials.
!-----------------------------------------------------------------------
IF(medium)THEN
   WRITE(Out_Unit,*)'Hermite polynomials evaluated at x = ', x
   WRITE(Out_Unit,*)'maximum order = ', norder
ENDIF
IF(full)THEN
   DO 5 n = 0, norder
      WRITE(Out_Unit,*)'n= ', n, '  Hn(x) = ', hn(n)
    5    CONTINUE
ENDIF



!
RETURN
END
