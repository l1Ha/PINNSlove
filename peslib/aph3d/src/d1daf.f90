REAL(Kind=WP_Kind) FUNCTION d1daf (x1, x2, nmx, sgmx)
USE Numeric_Kinds_Module
USE Numbers_Module
!* --------------------------------------------------------------
!* The DAF function of d^2/(dx^2) reads
!*
!*   d1daf(x-x')=-1/(2*sgmx^2*sqrt(pi))
!*             *sum_{n=0}^{M/2} (-1/4)^n/(n!)
!*             *H_(2n+2)((x-x')/(sqrt(2)*sgmx))
!*             *exp(-(x-x')^2/(2.0*sgmx^2))           .
!*
!* where H_(2n+2)( ) is a Hermite polynomial.
!*
!* On a grid, the kinetic energy DAF is
!*
!*   kdaf(i,j)=-dx/(2.0*m)*d1daf(x_i-x_j)
!*
!* where dx is the grid spacing and m is the mass.
!*
!* In the calculation, nmx=40, sgmx/dx = 2.85 are used.
!*
!* If you have any question, please let me know. My phone number
!* is 713-743-3247.
!* Youhong Huang
!* --------------------------------------------------------------
IMPLICIT NONE
INTEGER :: nmx, i, i1

REAL(Kind=WP_Kind) :: x1, x2, sgmx, sqt2, fac, w, dx, x, y, z, exp, hermite (0:1000)
sqt2 = SQRT (2.0d0)
fac = 1.0d0

w = 0.0d0

dx = x1 - x2
x = 1.0d0
hermite (0) = x
y = 2.0d0 * dx / (sqt2 * sgmx)

hermite (1) = y
DO i = 2, 2 * nmx + 2
i1 = i - 1
z = 2.0d0 * dx / (sqt2 * sgmx) * y - 2.0d0 * REAL(i1,WP_Kind) * x
hermite (i) = z
x = y
y = z
IF(mod (i, 2) ==0)THEN
   w = w + hermite (i - 1) * fac
   fac = fac * ( - 0.25d0) * 2.0d0 / REAL(i,WP_Kind)
ENDIF

ENDDO

d1daf = - w / (2.0d0 * sgmx**2.d0 * sqrt (pi) ) * exp ( - 0.5d0 * (dx / sgmx) **2.d0)
RETURN
END FUNCTION d1daf
