PROGRAM Program_TestReadAllData
USE FileUnits_Module
IMPLICIT NONE

OutDir=TRIM(OutputRoot)//"Testing/Test_ReadAllData/"
InputDir=TRIM(InputRoot)//"Testing/Test_ReadAllData/"

CALL ReadAllData

STOP 'TestReadAllData Completed Successfully'
END PROGRAM Program_TestReadAllData