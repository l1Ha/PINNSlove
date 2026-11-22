REAL(Kind=WP_Kind) FUNCTION rand(dummy)
!-----------------------------------------------------------------------
! PORT COMPAT (2026-10-01, gfortran port):
!   VMS/Alliant style random number generator used by RANFX.  Returns a
!   uniform deviate in [0,1).  The dummy argument is a control flag in the
!   original runtime; it is ignored (seeding is handled by the RNG state).
!-----------------------------------------------------------------------
USE Numeric_Kinds_Module
IMPLICIT NONE
INTEGER, INTENT(IN) :: dummy
REAL(Kind=WP_Kind) :: x

CALL RANDOM_NUMBER(x)
rand = x
RETURN
END FUNCTION rand

REAL(Kind=WP_Kind) FUNCTION secnds(x)
!-----------------------------------------------------------------------
! PORT COMPAT (2026-10-01, gfortran port):
!   VMS SECNDS: seconds since midnight; the argument is subtracted.  Only
!   used for timing printouts, so precision is unimportant.
!-----------------------------------------------------------------------
USE Numeric_Kinds_Module
IMPLICIT NONE
REAL(Kind=WP_Kind), INTENT(IN) :: x
INTEGER :: values(8)
REAL(Kind=WP_Kind) :: t

CALL DATE_AND_TIME(VALUES=values)
t = 3600.d0*values(5) + 60.d0*values(6) + values(7) + 1.d-3*values(8)
secnds = t - x
RETURN
END FUNCTION secnds
