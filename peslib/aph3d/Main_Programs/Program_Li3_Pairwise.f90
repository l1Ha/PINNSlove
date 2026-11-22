PROGRAM Program_Li3_Pairwise
USE FileUnits_Module
IMPLICIT NONE
CHARACTER(LEN=100), PARAMETER:: SysDir='Li3/APH3D_Li3_Pairwise/'

OutDir=TRIM(OutputRoot)//TRIM(SysDir)
InputDir=TRIM(InputRoot)//TRIM(SysDir)

CALL BasisForScattering

STOP 'Program_Li3_Pairwise'
ENDPROGRAM Program_Li3_Pairwise