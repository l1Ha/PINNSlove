PROGRAM Program_TestBasisRegions
USE FileUnits_Module
IMPLICIT NONE
LOGICAL CallFEM, CallDVR, CallPDAF, CallABM, CallDelves

CallFEm=.False.
CallPDAF=.False.
CallDelves=.False.
CallDVR=.True.
CallABM=.True.

OutDir=TRIM(OutputRoot)//"Testing/Test_BasisRegions/"
InputDir=TRIM(InputRoot)//"Testing/Test_BasisRegions/"

CALL ReadAllData
CALL BasisRegions(CallFEM, CallDVR, CallPDAF, CallABM, CallDelves)

STOP 'TestBasisRegions Completed Successfully'
END PROGRAM Program_TestBasisRegions