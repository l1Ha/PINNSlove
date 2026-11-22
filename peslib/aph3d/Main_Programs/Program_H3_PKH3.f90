PROGRAM Program_H3_PKH3
USE FileUnits_Common_Module
IMPLICIT NONE
Integer test
CHARACTER(LEN=100), PARAMETER:: SysDir='APH3D_H3_PKH3/'

OutDir=TRIM(OutputRoot)//"H3/"//TRIM(SysDir)
InputDir=TRIM(InputRoot)//"H3/"//TRIM(SysDir)
CALL BasisForScattering

STOP ' Program_H3_PKH3'
ENDPROGRAM Program_H3_PKH3