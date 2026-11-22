SUBROUTINE Laguer_fun ( y, n, h, m )
USE Numeric_Kinds_Module
USE FileUnits_Module

!
!--------------------- Laguer ------------------------------------------
!
! Laguer -- computes Laguere polynomials
!
! --- on input the user supplies --
!
!     y  -- the dimensionless argument of the Laguere polynomial
!     n  -- the order of the polynomial (n=0 ==> h0=1)
!     m  -- n+1 (required), also the length of the space in h
!
! --- on exit the routine returns --
!
!     h  -- vector of length m, containing h(1)=h0 h(2)=h1, etc.
!     Laguer -- function value returned is h(m) = nth Laguere poly.
!
!-----------------------------------------------------------------------
!
INTEGER n, m, i
REAL(Kind=WP_Kind) h(0:n), y
!
!-------------------- execution begins here ----------------------------
!
h(0) = 1d0
IF( n == 0 ) RETURN
!
h(1) = 1.d0 - y
IF( n == 1 ) RETURN
!
DO 1 i  = 1, n-1
   h(i+1) = ((2*i+1-y)*h(i) - i*h(i-1))/(i+1)
 1    CONTINUE

RETURN
END
