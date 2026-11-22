PROGRAM Program_HydrogenPositron
USE FileUnits_Module
IMPLICIT NONE
LOGICAL, PARAMETER:: H3_Calc=.True.
CHARACTER(LEN=100), PARAMETER:: SysDir='APH3D_H_Positron/'

OutDir=TRIM(OutputRoot)//'H+Particle/'//TRIM(SysDir)
InputDir=TRIM(InputRoot)//'H+Particle/'//TRIM(SysDir)


CALL BasisForScattering

STOP 'Program_HydrogenPositron'
ENDPROGRAM Program_HydrogenPositron