SUBROUTINE alagur(lk, kmax, z, alag)
USE Numeric_Kinds_Module
!

!  calculates the lowest kmax associated laguerre polynomials by using
!  recurrence relation (13.40) from arfken, mathematical methods for
!  physicists, academic press, 1970.  to use for an hydrogen atom,
!  lk = 2l+1, k = n - l - 1, z = ar where a = 2/n a.u.
!-----------------------------------------------------------------------
!  written by b.j. archer
!-----------------------------------------------------------------------
USE FileUnits_Module
IMPLICIT NONE
REAL(Kind=WP_Kind) alag, z

INTEGER k, kmax, lk

INTRINSIC dfloat

DIMENSION alag(0:100)

DO k=0, 100
   alag(k)=0.0d0
ENDDO

alag(0) = 1.0d0
alag(1) = lk + 1.0d0 - z
DO k=2, kmax
   alag(k) = ((2.0d0*k+lk-1.0d0-z)*alag(k-1) -(k+lk-1.0d0)*alag(k-2))/(1.d0*k)
ENDDO

RETURN
!---------------***END-alagur***----------------------------------------
ENDSUBROUTINE alagur
