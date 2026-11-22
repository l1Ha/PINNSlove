PROGRAM Program_HydrogenElectron
USE FileUnits_Module
IMPLICIT NONE
LOGICAL, PARAMETER:: H3_Calc=.True.
CHARACTER(LEN=100), PARAMETER:: SysDir='APH3D_H_Electron\'

OutDir=TRIM(OutputRoot)//'H+Particle/'//TRIM(SysDir)
InputDir=TRIM(InputRoot)//'H+Particle/'//TRIM(SysDir)


CALL BasisForScattering

STOP 'Program_HydrogenElectron'
ENDPROGRAM Program_HydrogenElectron