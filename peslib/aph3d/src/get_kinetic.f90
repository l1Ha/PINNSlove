SUBROUTINE get_kinetic (n, toe, kinetic, xvals, kmax, kase, system, group, tblock, iblock, usys2, dafx, nd_daf, ndaf)
USE Numeric_Kinds_Module
USE FileUnits_Module
IMPLICIT NONE
LOGICAL, PARAMETER :: debug=.False.
INTEGER :: i, j, n, kmax
INTEGER :: nd_daf, ndaf, iblock (n, n)
CHARACTER(LEN=3) :: system, group
CHARACTER(LEN=8) :: kase
REAL(Kind=WP_Kind) :: toe (0:kmax), kinetic (n, n), xvals (n), tblock (n, n), usys2, dafx (0:nd_daf, 0:ndaf)
WRITE(Out_Unit, * )
WRITE(Out_Unit, * ) 'Begin(Get_Kinetic)', kase

DO i = 1, n
   DO j = 1, n
      kinetic (i, j) = 0.d0
   ENDDO
   IF(kase=='chi     ') tblock (1, 1)  = 1.d0
   IF(kase=='coriolis') tblock (1, 1)  = 1.d0
ENDDO

IF(n>1)THEN
   IF(kase=='rho     ')THEN
      CALL kin_rho (n, toe, kinetic, xvals, kmax)
   ELSEIF (kase=='theta   ')THEN
      CALL kin_theta (n, kinetic, xvals, kmax, usys2, dafx, nd_daf, ndaf)
   ELSEIF (kase=='chi     ')THEN
      CALL kin_chi (n, toe, kinetic, xvals, kmax, system, group, tblock, iblock)
   ELSEIF (kase=='coriolis')THEN
      WRITE(Out_Unit, * ) 'coriolis'
      CALL kin_chi (n, toe, kinetic, xvals, kmax, system, group, tblock, iblock)
   ELSE
      WRITE(Out_Unit, * ) 'Incorrect kase=', kase
      STOP 'get_kinetic'
   ENDIF
ELSE
   kinetic(1,1)=1.d0
ENDIF

IF(debug)THEN
   WRITE(Out_Unit, * ) 'Kinetic Energy Matrix, kase=', kase
   CALL MxOut(kinetic, n, n)
ENDIF

WRITE(Out_Unit, * ) 'END(Get_Kinetic)', kase

WRITE(Out_Unit, * )
RETURN
ENDSUBROUTINE get_kinetic
