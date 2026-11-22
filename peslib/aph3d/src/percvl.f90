REAL(Kind=WP_Kind) FUNCTION percvl (lam, j, l, jp, lp, jtot)
!
!  REAL(Kind=WP_Kind) function percvl (lam, j, l, jp, lp, jtot)
USE FileUnits_Module
IMPLICIT NONE
INTEGER lam, j, l, jp, lp, jtot
REAL(Kind=WP_Kind) c1, c2, cleb1, cleb2, factor, racah, c6j
!#include <Inc.unitno>
!-----------------------------------------------------------------------
!  calculates the percival seaton coefficients.
!-----------------------------------------------------------------------
c1=(-1.d0)**jp*sqrt((2*lam+1.d0)/(2*j+1.d0))*cleb1(lam, jp, j, 0, 0, 0)
c2=(-1.d0)**lp*sqrt((2*lam+1.d0)/(2*l+1.d0))*cleb2(lam, lp, l, 0, 0, 0)
factor=(-1)**(j+jp-jtot)/(2*lam+1.d0)*sqrt((2*j+1.d0)*(2*l+1.d0)*(2*jp+1.d0)*(2*lp+1.d0))
racah=c6j(lam, l, lp, jtot, jp, j)*(-1)**(j+l+jp+lp)
percvl=factor*c1*c2*racah
IF(ABS(percvl) < 1.d-12) percvl=0.d0
RETURN
END
