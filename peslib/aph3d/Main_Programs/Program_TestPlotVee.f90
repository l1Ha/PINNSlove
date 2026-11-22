PROGRAM Program_Test_PlotVee
USE Numeric_Kinds_Module
USE FileUnits_Common_Module
USE InputFile_Module
USE BasisFlowOptions_Module
IMPLICIT NONE
LOGICAL, PARAMETER:: H3_Calc=.True.
CHARACTER(LEN=10) hour
CHARACTER(LEN=8) today
CHARACTER(LEN=5) curzone

OUtDir=TRIM(OutputRoot)//"Testing/Test_PlotVee/"
InputDir=TRIM(InputRoot)//"Testing/Test_PlotVee/"

!OPEN(Unit=DBug_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'Output/debuging.txt',Form='FORMATTED')
OPEN(Unit=Out_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'Output/BasisForScatt.txt',Form='FORMATTED')

CALL Date_And_Time(today, hour, curzone, dtvalues)                  ! Write Time and Data
WRITE(Out_Unit,'("At the start of AAA_Basis_Main")')
CALL Print_Date_And_Time(Today, Hour, CurZone, DtValues, Out_Unit)

WRITE(Out_Unit,*)'Read in All of the Input Data'
CALL ReadAllData
WRITE(Out_Unit,*)'Determine Basis Set Distances'
CALL BasisRegions(CallFEM, CallDVR, CallPDAF, CallABM, CallDelves)

PlotVee=.True.

IF(PlotVee)THEN
   WRITE(Out_Unit,*)'Plot Colinear Potential'
   WRITE(Msg_Unit,*)'Plot Colinear Potential'
   CALL Graf_Collinear(Basis_Dist(NSectors))
   
   WRITE(Out_Unit,*)'Plot Planar Slices of the Potential'
   WRITE(Msg_Unit,*)'Plot Planar Slices of the Potential'
   CALL Graf_Planar(Basis_Dist(NSectors))
   
   WRITE(Out_Unit,*)'Plot Hypersurfaces of the Potential'
   WRITE(Msg_Unit,*)'Plot Hypersurfaces of the Potential'
   Call Graf_Hyper
   
   WRITE(Out_Unit,*)'Plot Potential in Triangular Coordinates'
   WRITE(Msg_Unit,*)'Plot Potential in Triangular Coordinates'
   CALL Graf_Triangular
   
   WRITE(Out_Unit,*)'Plot Potential Minima'
   WRITE(Msg_Unit,*)'Plot Potential Minima'
   CALL Graf_Surm
ENDIF

STOP 'Test_PlotVee'
END PROGRAM Program_Test_PlotVee