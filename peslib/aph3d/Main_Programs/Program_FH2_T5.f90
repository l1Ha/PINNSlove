PROGRAM Program_FH2_T5
USE FileUnits_Module
IMPLICIT NONE
CHARACTER(LEN=100), PARAMETER:: SysDir='APH3D_FH2_T5/'

OutDir=TRIM(OutputRoot)//'FH2/'//TRIM(SysDir)
InputDir=TRIM(InputRoot)//'FH2/'//TRIM(SysDir)

CALL BasisForScattering

STOP 'Program_FH2_T5'
END PROGRAM Program_FH2_T5