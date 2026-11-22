SUBROUTINE symham (n, hamil, tblock, iblock, system, group)
USE Numeric_Kinds_Module
USE FileUnits_Module
IMPLICIT NONE
LOGICAL :: debug
CHARACTER(LEN=3) :: group, system
INTEGER :: n, nirrep, gamma, irrep, j, i, nsymop, iw, isymop
INTEGER :: iblock (n, n), isum, k, n1 (8), nlast, index, fac, par (8)
INTEGER :: nb, ne

REAL(Kind=WP_Kind) :: hamil (n, n), tblock (n, n), norm (1000), cpu_begin, cpu_end
debug = .true.
WRITE(Out_Unit, * ) 'system=', system

WRITE(Out_Unit, * ) 'group=', group
CALL cputime (cpu_begin)
DO i = 1, n
   DO j = 1, n
      iblock (i, j) = 0
      tblock (i, j) = 0.d0
   ENDDO
ENDDO
CALL cputime (cpu_end)

WRITE(Out_Unit, * ) 'block1: time=', cpu_end-cpu_begin
nlast = 0
CALL symop (irrep, 1, i, group, gamma, n, isymop, n1, nirrep, nsymop, 0, system, fac, par)
DO irrep = 1, nirrep
   CALL cputime (cpu_begin)
   iw = 0
   nb = nlast + 1
   ne = nlast + n1 (irrep)
   DO j = nb, ne
      isum = 0
      DO WHILE (isum==0)
         DO isymop = 1, nsymop
            index = j - nlast + iw
CALL symop (irrep, index, i, group, gamma, n, isymop, n1, nirrep, nsymop, 1, system, fac, par)
            iblock (i, j) = gamma + iblock (i, j)
            IF(fac/=0.AND.i/=1)THEN
               iblock (n + 2 - i, j) = iblock (n + 2 - i, j) + fac * gamma
            ELSEIF (i==1)THEN
               iblock (i, j) = iblock (i, j) + fac * gamma
            ENDIF
         ENDDO
         isum = 0
         DO k = 1, n
            isum = isum + iblock (k, j) * iblock (k, j)
         ENDDO
         IF(isum==0) iw = iw + 1
      ENDDO
      IF(j/=nlast + 1)THEN
   CALL ischmidt (iblock (1, nlast + 1), n, j - nlast - 1, iblock (1, j) )
      ENDIF
   ENDDO
CALL cputime (cpu_end)

WRITE(Out_Unit, * ) 'block2: time=', cpu_end-cpu_begin
CALL cputime (cpu_begin)
nlast = nlast + n1 (irrep)
DO i = nb, ne
DO j = i, ne
isum = 0
DO k = 1, n
isum = isum + iblock (k, i) * iblock (k, j)
ENDDO
ENDDO
ENDDO
CALL cpu_time (cpu_end)
WRITE(Out_Unit, * ) 'block3: time=', cpu_end-cpu_begin

ENDDO
CALL cputime (cpu_begin)
DO i = 1, n
!ccc         DO j=i,n
j = i
isum = 0
DO k = 1, n
isum = isum + iblock (k, i) * iblock (k, j)
ENDDO
IF(i==j)THEN
   IF(isum==0)THEN
      WRITE(Out_Unit, * ) 'error ,i,j,isum=', i, j, isum
   ENDIF
   norm (i) = isum
   norm (i) = 1.d0 / sqrt (norm (i) )
ELSE
   IF(isum/=0)THEN
      WRITE(Out_Unit, * ) 'error', i, j, isum
   ENDIF
ENDIF
!ccc         ENDDO

ENDDO


IF(debug)THEN
   WRITE(Out_Unit, * ) 'iblock matrix'
   DO i = 1, min (n, 25)
   WRITE(Out_Unit, 10) (iblock (i, j), j = 1, min (n, 25) )
   ENDDO
ENDIF
   10 FORMAT(1x,25i3)
CALL cputime (cpu_end)

WRITE(Out_Unit, * ) 'block4: time=', cpu_end-cpu_begin
CALL cputime (cpu_begin)
DO i = 1, n
   DO j = 1, n
      tblock (i, j) = 0.d0
      DO k = 1, n
         tblock (i, j) = tblock (i, j) + hamil (i, k) * iblock (k, j)
      ENDDO
   ENDDO
ENDDO
DO i = 1, n
   DO j = 1, n
      hamil (i, j) = tblock (i, j) * norm (j)
   ENDDO
ENDDO
DO i = 1, n
   DO j = 1, n
      tblock (i, j) = 0.d0
      DO k = 1, n
         tblock (i, j) = tblock (i, j) + iblock (k, i) * hamil (k, j)
      ENDDO
   ENDDO
ENDDO
DO i = 1, n
   DO j = 1, n
      hamil (i, j) = tblock (i, j) * norm (i)
   ENDDO
ENDDO
CALL cputime (cpu_end)

WRITE(Out_Unit, * ) 'block5: time=', cpu_end-cpu_begin
CALL cputime (cpu_begin)
nlast = 1
DO irrep = 1, nirrep
   DO i = 1, n
      DO j = 1, n
         IF(i>=nlast.AND.i<nlast + n1 (irrep) )THEN
            IF(j<nlast.or.j>=nlast + n1 (irrep) )THEN
               IF(abs (hamil (i, j) ) <1.d-10)THEN
                  hamil (i, j) = 0.d0
               ELSE
         WRITE(Out_Unit, * ) 'irrep,i,j=', irrep, i, j
         WRITE(Out_Unit, * ) 'Error hamil=', hamil (i, j),ABS(hamil (i, j)),1.d-10,abs (hamil (i, j) ) <1.d-10
                  STOP 'symham'
               ENDIF
            ENDIF
         ENDIF
   ENDDO
ENDDO
nlast = nlast + n1 (irrep)

ENDDO
IF(debug)THEN
   WRITE(Out_Unit, * ) 'Hamiltonian matrix'
   CALL MxOut(hamil, n, n)
ENDIF
   20 FORMAT(1x,24f5.0)
CALL cputime (cpu_end)

WRITE(Out_Unit, * ) 'block6: time=', cpu_end-cpu_begin
DO i = 1, n
   DO j = 1, n
      tblock (i, j) = iblock (i, j) * norm (j)
   ENDDO
ENDDO

RETURN
ENDSUBROUTINE symham
