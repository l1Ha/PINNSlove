PROGRAM Program_H3_DMBE
USE FileUnits_Module
LOGICAL, PARAMETER:: H3_Calc=.True.
CHARACTER(LEN=100), PARAMETER:: SysDir='APH3D_H3_DMBE/'

CALL APH3D_InitRoots

OutDir=TRIM(OutputRoot)//TRIM(SysDir)
InputDir=TRIM(InputRoot)//TRIM(SysDir)

CALL BasisForScattering

STOP 'Program_H3_DMBE'
END PROGRAM Program_H3_DMBE