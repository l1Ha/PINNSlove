PROGRAM Program_Li2
USE FileUnits_Module
IMPLICIT NONE
CHARACTER(LEN=100), PARAMETER:: SysDir='APH3D_Li2/'

OutDir=TRIM(OutputRoot)//TRIM(SysDir)
InputDir=TRIM(InputRoot)//TRIM(SysDir)

CALL BasisForScattering

STOP 'Program_Li2'
ENDPROGRAM Program_Li2