PROGRAM Program_He3_CC
USE FileUnits_Module
IMPLICIT NONE
LOGICAL, PARAMETER:: H3_Calc=.True.
CHARACTER(LEN=100), PARAMETER:: SysDir='He3/APH3D_He3_CC/'

OutDir=TRIM(OutputRoot)//TRIM(SysDir)
InputDir=TRIM(InputRoot)//TRIM(SysDir)

CALL BasisForScattering

STOP 'Program_He3_CC'
END PROGRAM Program_He3_CC