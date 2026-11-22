PROGRAM Program_LiH2
USE FileUnits_Module
IMPLICIT NONE
CHARACTER(LEN=100), PARAMETER:: SysDir='APH3D_LiH2/'

OutDir=TRIM(OutputRoot)//"LiH2/"//TRIM(SysDir)
InputDir=TRIM(InputRoot)//"LiH2/"//TRIM(SysDir)

CALL BasisForScattering

STOP 'Program_LiH2'
ENDPROGRAM Program_LiH2