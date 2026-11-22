SUBROUTINE Hermit_Fun ( y, n, h, m )
!
!--------------------- Hermit ------------------------------------------
!
! Hermit -- computes Hermite polynomials
!
! --- on input the user supplies --
!
!     y  -- the dimensionless argument of the Hermite polynomial
!     n  -- the order of the polynomial (n=0 ==> h0=1)
!     m  -- n+1 (required), also the length of the space in h
!
! --- on exit the routine returns --
!
!     h  -- vector of length m, containing h(1)=h0 h(2)=h1, etc.
!     Hermit -- function value returned is h(m) = nth Hermite poly.
!
!----------------------------------------------------------------------
!
USE Numeric_Kinds_Module
USE FileUnits_Module
IMPLICIT NONE
INTEGER n, m, ip, i
REAL(Kind=WP_Kind) h(m), y, ty, tn, Hermit
!
!------------------- execution begins here ----------------------------
!
h(1) = 1d0
IF( n == 0 ) RETURN
!
ty = 2d0*y
h(2) = ty
IF( n == 1 ) RETURN
!
tn = 0d0
!
DO ip  = 3, m
   i     = ip-1
   tn    = tn + 2d0
   h(ip) = ty*h(i) - tn*h(i-1)
ENDDO
!
Hermit = h(m)
!
RETURN
ENDSUBROUTINE Hermit_Fun
