SUBROUTINE Pn ( y, n, p)
!
!----------------------- Pn ------------------------------------------
!
! Pn -- computes Legendre  polynomials
!
! --- on input the user supplies --
!
!     y  -- the dimensionless argument of the Legendre  polynomial
!     n  -- the order of the polynomial (n=0 ==> p0=1)
!
! --- on exit the routine returns --
!
!     p  -- vector of length (n+1), containing p(0)=p0 p(1)=p1, etc.
!     Legendre -- function value returned is p(n) = nth Legendre  poly.
!
!----------------------------------------------------------------------
!
USE Numeric_Kinds_Module
USE FileUnits_Module
IMPLICIT NONE
INTEGER n, i
REAL(Kind=WP_Kind) y, p(0:n)
!
!------------------- execution begins here ----------------------------
!
p(0) = 1d0
IF( n == 0 ) RETURN
!
p(1) = y
IF( n == 1 ) RETURN
!
DO  i  = 1, n-1
   p(i+1) = y*p(i) + (i*(y*p(i)-p(i-1)))/(i+1)
ENDDO
!
!   Pn = p(n)
!
RETURN
ENDSUBROUTINE Pn
