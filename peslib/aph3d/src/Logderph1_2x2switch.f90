SUBROUTINE Get2x2Switch(ps2x2)
!-----------------------------------------------------------------------
! PORT (2026-10-02): reads the APH3D_PS2X2 environment switch used by the
! two-electronic-state (2x2) data path of the inner-region propagation.
!-----------------------------------------------------------------------
USE Numeric_Kinds_Module
IMPLICIT NONE
LOGICAL, INTENT(OUT) :: ps2x2
CHARACTER(LEN=64) :: env2
INTEGER :: len2, stat2
ps2x2 = .FALSE.
CALL get_environment_variable('APH3D_PS2X2', env2, len2, stat2)
IF(stat2==0 .AND. len2>0)THEN
   IF(env2(1:1)=='1' .OR. env2(1:1)=='T' .OR. env2(1:1)=='t') ps2x2=.TRUE.
ENDIF
RETURN
ENDSUBROUTINE Get2x2Switch
