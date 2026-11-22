PROGRAM Program_HeH2
USE FileUnits_Module
IMPLICIT NONE
LOGICAL, PARAMETER:: H3_Calc=.True.
CHARACTER(LEN=100), PARAMETER:: SysDir='HeH2/APH3D_HeH2/'

OutDir=TRIM(OutputRoot)//TRIM(SysDir)
InputDir=TRIM(InputRoot)//TRIM(SysDir)

CALL BasisForScattering

STOP 'Program_HeH2'
ENDPROGRAM Program_HeH2