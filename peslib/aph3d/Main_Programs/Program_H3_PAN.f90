PROGRAM Program_H3_PAN
USE FileUnits_Module
IMPLICIT NONE
LOGICAL, PARAMETER:: H3_Calc=.True.
CHARACTER(LEN=100), PARAMETER:: SysDir='APH3D_H3_PAN/'

OutDir=TRIM(OutputRoot)//TRIM(SysDir)
InputDir=TRIM(InputRoot)//TRIM(SysDir)


CALL BasisForScattering

STOP 'Program_H3_PAN'
ENDPROGRAM Program_H3_PAN