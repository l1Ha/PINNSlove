PROGRAM Program_TestDelvesiBasis
USE FileUnits_Common_Module
IMPLICIT NONE

OutDir=TRIM(OutputRoot)//"Testing/Test_DelvesBasis/"
InputDir=TRIM(InputRoot)//"Testing/Test_DelvesBasis/"

CALL ReadAllData
CALL Delves_Basis

STOP 'TestDelvesiBasis Completed Successfully'
END PROGRAM Program_TestDelvesiBasis