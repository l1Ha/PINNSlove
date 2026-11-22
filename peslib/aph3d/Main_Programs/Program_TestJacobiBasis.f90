PROGRAM Program_TestJacobiBasis
USE FileUnits_Common_Module
IMPLICIT NONE

OutDir=TRIM(OutputRoot)//"Testing/Test_JacobiBasis/"
InputDir=TRIM(InputRoot)//"Testing/Test_JacobiBasis/"

CALL ReadAllData
CALL Jacobi_Basis

STOP 'TestJacobiBasis Completed Successfully'
END PROGRAM Program_TestJacobiBasis