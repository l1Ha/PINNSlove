SUBROUTINE hermit_ex (y, n, h, expfac)
USE Numeric_Kinds_Module
USE Numbers_Module
!-----------------------------------------------------------------------
! This routine was written by G. A. Parker
! If you find an error or have an improvement please send a messge to
! parker@ou.edu
!-----------------------------------------------------------------------

!--------------------- hermit ------------------------------------------
!
! >>> hermit -- computes hermite polynomials
!
! --- on input the user supplies --
!
!     y  -- the dimensionless argument of the hermite polynomial
!     n  -- the order of the polynomial (n=0 ==> h0=1)
!
! --- on exit the routine returns --
!
!     h  -- vector of length m, containing h(0)=h0 h(1)=h1, etc.
!
!----------------------------------------------------------------------

IMPLICIT NONE

INTEGER i, n
REAL(KIND=WP_Kind) expfac, y, ty, tn
REAL(KIND=WP_Kind) h(0:n+1)

!------------------- execution begins here ----------------------------

      h(0) = One*expfac
      IF( n .eq. 0 ) RETURN

      ty = Two*y
      h(1) = ty*expfac
      IF( n .eq. 1 ) RETURN

DO i  = 1, n-1, 2
   h(i+1) = -y*h(i)/(i+1) + i*h(i-1)/(i+1)
   h(i+2) =  ty*h(i+1) + h(i)
ENDDO

RETURN
ENDSUBROUTINE hermit_ex
