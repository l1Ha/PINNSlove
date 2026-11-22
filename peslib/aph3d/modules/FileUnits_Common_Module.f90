MODULE FileUnits_Common_Module
   USE Numeric_Kinds_Module
   IMPLICIT NONE
   SAVE
   CHARACTER(LEN=99)  OutDir
   CHARACTER(LEN=99)  InputDir
   ! PORT NOTE (2026-10-01): roots are runtime configurable so the same build
   ! runs on any machine.  Call APH3D_InitRoots (from the main program, before
   ! OutDir/InputDir are set) or export APH3D_INPUT_ROOT / APH3D_OUTPUT_ROOT.
   CHARACTER(LEN=99) :: InputRoot = "D:/Input_Data/APH3D/"
   CHARACTER(LEN=99) :: OutputRoot = "E:/ParkerE/APH3D/"
   !INTEGER(KIND=IW_Kind), PARAMETER::Screen_Unit   =  0
   INTEGER(KIND=IW_Kind), PARAMETER::Msg_Unit      =  0        !
   INTEGER(KIND=IW_Kind), PARAMETER::In_Unit       =  3        !
   INTEGER(KIND=IW_Kind), PARAMETER::Out_Unit      =  6        !
CONTAINS
   SUBROUTINE APH3D_InitRoots()
      CHARACTER(LEN=256):: envval
      INTEGER:: leng, stat
      CALL GET_ENVIRONMENT_VARIABLE('APH3D_INPUT_ROOT', envval, leng, STATUS=stat)
      IF(stat==0 .AND. leng>0) InputRoot=envval(1:leng)
      CALL GET_ENVIRONMENT_VARIABLE('APH3D_OUTPUT_ROOT', envval, leng, STATUS=stat)
      IF(stat==0 .AND. leng>0) OutputRoot=envval(1:leng)
   ENDSUBROUTINE APH3D_InitRoots
ENDMODULE FileUnits_Common_Module
