SUBROUTINE cputime (cpu)
USE Numeric_Kinds_Module

IMPLICIT NONE
REAL(Kind=WP_Kind) :: cpu
real :: tclock (2)
EXTERNAL etime

!c      CALL etime(tclock)
cpu = tclock (1)

cpu = 0.d0
RETURN
ENDSUBROUTINE cputime
