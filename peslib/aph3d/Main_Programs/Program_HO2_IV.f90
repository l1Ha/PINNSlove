PROGRAM Program_HO2_IV
USE FileUnits_Module
IMPLICIT NONE
LOGICAL, PARAMETER:: H3_Calc=.True.
CHARACTER(LEN=100), PARAMETER:: SysDir='HO2/APH3D_HO2_IV/'

OutDir=TRIM(OutputRoot)//TRIM(SysDir)
InputDir=TRIM(InputRoot)//TRIM(SysDir)


CALL BasisForScattering

STOP 'Program_HO2_IV'
END PROGRAM Program_HO2_IV