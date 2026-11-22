PROGRAM Program_Rb2
USE FileUnits_Module
IMPLICIT NONE
CHARACTER(LEN=100), PARAMETER:: SysDir='APH3D_Rb2/'

OutDir=TRIM(OutputRoot)//TRIM(SysDir)
InputDir=TRIM(InputRoot)//TRIM(SysDir)

CALL BasisForScattering

STOP 'Program_Rb2'
ENDPROGRAM Program_Rb2