SUBROUTINE interpfn (h, p, xpt, wht, nquad)
USE Numeric_Kinds_Module
USE Numbers_Module
!
!   $RCSfile: interpfn.f,v $ $Revision: 1.15 $
!   $Date: 89/11/16 11:02:59 $
!   $State: Stable $
!

!              P U R P O S E   O F   S U B R O U T I N E

!    Supplies the coefficients for evaluating the surface functions
!    and their derivitives.

!                  I N P U T   A R G U M E N T S

!  xpt
!  wht
!  nquad


!               O U T P U T   A R G U M E N T S

!    h        Interpolation coefficients for evaluating the surface
!             functions.
!    p        Interpolation coefficients for evaluating the derivatives
!             of the surface functions.

IMPLICIT NONE
INTEGER i, j, nquad
REAL(Kind=WP_Kind) r, s, xpt, wht, h, p, xj, det, xx
DIMENSION h(9, nquad, nquad), p(2,9, nquad, nquad), xj(2,2),  xx(2,9), xpt(nquad), wht(nquad)
EXTERNAL funct2, glegen
!-----------------------------------------------------------------------
! Get Gauss_Legendre quadrature points and abscissa
!-----------------------------------------------------------------------
CALL glegen (nquad, xpt, wht, -One, One)

!-----------------------------------------------------------------------
! calculate interpolation functions.
!-----------------------------------------------------------------------
DO 25 i = 1, nquad
   r = xpt(i)
   DO 20 j = 1, nquad
      s = xpt(j)
      CALL funct2 (r, s, h(1, i, j), p(1, 1, i, j), xj, det, xx, 1)
   20    CONTINUE
   25 CONTINUE

RETURN
!
END
