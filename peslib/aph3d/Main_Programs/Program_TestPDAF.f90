PROGRAM Program_Test_PDAF
USE FileUnits_Common_Module
IMPLICIT NONE
LOGICAL, PARAMETER:: H3_Calc=.True.

OutDir="E:/ParkerE/APH3D/Test_PDAF/"
InputDir="C:/Users/Parker/Documents/Input_Data/APH3D/Test_PDAF/"

CALL PDAF

STOP 'Test_Asymptotic'
END PROGRAM Program_Test_PDAF