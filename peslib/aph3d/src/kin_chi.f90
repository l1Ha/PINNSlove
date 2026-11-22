SUBROUTINE kin_chi (n, toe, kinetic, xvals, kmax, system, group, tblock, iblock)
USE Numeric_Kinds_Module
USE FileUnits_Module
IMPLICIT NONE
LOGICAL :: debug
CHARACTER(LEN=3) :: system, group
INTEGER :: i, j, n, m, index, kmax, mloop, iblock (n, n)
REAL(Kind=WP_Kind) :: toe (0:kmax), kinetic (n, n), xvals (n), tblock (n, n), cpu_begin, cpu_end
WRITE(Out_Unit, * )
WRITE(Out_Unit, * ) 'Begin(Kin_Chi)'

debug = .true.
mloop = kmax / n + 1
CALL cputime (cpu_begin)
IF(debug) WRITE(Out_Unit, * ) 'mloop=', mloop, n, kmax
DO m = - mloop, mloop
   DO i = 1, n
      DO j = 1, n
         index = iabs (i - j + m * n)
         IF(index<kmax)THEN
            kinetic (i, j) = kinetic (i, j) + toe (index)
         ENDIF
      ENDDO
   ENDDO
ENDDO

CALL cputime (cpu_end)

WRITE(Out_Unit, * ) 'mloop: time=', cpu_end-cpu_begin
IF(debug)THEN
   WRITE(Out_Unit, * ) 'Chi-Kinetic Energy (before trans)'
   CALL MxOut(kinetic, n, n)
ENDIF


DO i=1,n
   DO j=1,i
      IF(ABS(kinetic(i,j)-kinetic(j,i))>1.0d-8)THEN
         WRITE(Msg_Unit,*)'WARNING: Hamil_chi not sym',i,j,kinetic(i,j),kinetic(j,i)
         STOP 'wrong in hamil_chi'
      ENDIF
   ENDDO
ENDDO


CALL cputime (cpu_begin)
CALL symham (n, kinetic, tblock, iblock, system, group)
CALL cputime (cpu_end)

WRITE(Out_Unit, * ) 'symham: time=', cpu_end-cpu_begin
IF(debug)THEN
   WRITE(Out_Unit, * ) 'Chi-Kinetic Energy (after trans)'
   CALL MxOut(kinetic, n, n)
ENDIF
WRITE(Out_Unit, * ) 'END(Kin_Chi)'

WRITE(Out_Unit, * )
RETURN
ENDSUBROUTINE kin_chi
