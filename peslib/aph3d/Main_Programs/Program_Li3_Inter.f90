PROGRAM Program_Li3_Inter
USE FileUnits_Module
IMPLICIT NONE
CHARACTER(LEN=100), PARAMETER:: SysDir='APH3D_Li3_Inter/'

OutDir=TRIM(OutputRoot)//TRIM(SysDir)
InputDir=TRIM(InputRoot)//TRIM(SysDir)

CALL BasisForScattering

STOP 'Program_Li3_Inter'
ENDPROGRAM Program_Li3_Inter