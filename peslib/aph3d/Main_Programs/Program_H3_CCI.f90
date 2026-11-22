PROGRAM Program_H3_CCI
USE FileUnits_Module
IMPLICIT NONE
LOGICAL, PARAMETER:: H3_Calc=.True.
CHARACTER(LEN=100), PARAMETER:: SysDir='APH3D_H3_CCI/'

OutDir=TRIM(OutputRoot)//TRIM(SysDir)
InputDir=TRIM(InputRoot)//TRIM(SysDir)

CALL BasisForScattering

STOP 'Program_H3_CCI'
END PROGRAM Program_H3_CCI