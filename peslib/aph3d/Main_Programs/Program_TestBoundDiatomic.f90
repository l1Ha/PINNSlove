PROGRAM Program_Test_BoundDiatomic
USE Numeric_Kinds_Module
USE FileUnits_Module
USE BasisFlowOptions_Module
IMPLICIT NONE
LOGICAL, PARAMETER:: H3_Calc=.True.
CHARACTER(LEN=10) hour
CHARACTER(LEN=8) today
CHARACTER(LEN=5) curzone

OutDir=TRIM(OutputRoot)//"Testing/Test_BoundDiatomic/"
InputDir=TRIM(InputRoot)//"Testing/Test_BoundDiatomic/"

OPEN(Unit=DBug_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'Output/debuging.txt',Form='FORMATTED')
OPEN(Unit=Out_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'Output/BasisForScatt.txt',Form='FORMATTED')

CALL Date_And_Time(today, hour, curzone, dtvalues)                  ! Write Time and Data
WRITE(Out_Unit,'("At the start of AAA_Basis_Main")')
CALL Print_Date_And_Time(Today, Hour, CurZone, DtValues, Out_Unit)

WRITE(Out_Unit,*)'Read in All of the Input Data'
CALL ReadAllData
CALL BoundDiatomic

STOP 'Test_BoundDiatomic'
END PROGRAM Program_Test_BoundDiatomic