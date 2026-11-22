PROGRAM Program_LiFH
USE FileUnits_Module
IMPLICIT NONE
CHARACTER(LEN=100), PARAMETER:: SysDir='APH3D_LiFH/'

OutDir=TRIM(OutputRoot)//"LiFH/"//TRIM(SysDir)
InputDir=TRIM(InputRoot)//"LiFH/"//TRIM(SysDir)

CALL BasisForScattering

STOP 'Program_LiFH'
ENDPROGRAM Program_LiFH