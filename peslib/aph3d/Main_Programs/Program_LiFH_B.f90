PROGRAM Program_LiFH_B
USE FileUnits_Module
IMPLICIT NONE
CHARACTER(LEN=100), PARAMETER:: SysDir='APH3D_LiFH_B/'
INTEGER jth, kth
REAL*8 E_nu_j, Yjk(4,4), term

OutDir=TRIM(OutputRoot)//"LiFH/"//TRIM(SysDir)
InputDir=TRIM(InputRoot)//"LiFH/"//TRIM(SysDir)

CALL BasisForScattering

STOP 'Program_LiFH_B'
ENDPROGRAM Program_LiFH_B