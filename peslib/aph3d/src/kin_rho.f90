SUBROUTINE kin_rho (n, toe, kinetic, xvals, kmax)
USE Numeric_Kinds_Module
USE FileUnits_Module
IMPLICIT NONE
LOGICAL :: debug
INTEGER :: i, j, n, index, kmax


REAL(Kind=WP_Kind) :: toe (0:kmax), kinetic (n, n), xvals (n), cpu_begin, cpu_end
WRITE(Out_Unit, * )

WRITE(Out_Unit, * ) 'Begin(Kin_Rho)'

debug = .False.
CALL cputime (cpu_begin)
DO i = 1, n
   DO j = 1, n
      index = iabs (i - j)
      IF(index<kmax)THEN
         kinetic (i, j) = kinetic (i, j) + toe (index)
      ENDIF
   ENDDO
ENDDO

CALL cputime (cpu_end)
IF(debug)THEN
   WRITE(Out_Unit, * ) 'Rho-Kinetic Energy Matrix'
   CALL MxOut(kinetic, n, n)
ENDIF
WRITE(Out_Unit, * ) 'END(Kin_Rho)'

WRITE(Out_Unit, * )
RETURN
ENDSUBROUTINE kin_rho
