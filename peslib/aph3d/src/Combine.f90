SUBROUTINE Combine
USE FileUnits_Module
USE BasisFlowOptions_Module, ONLY: CallFEM,CallDVR,CallPDAF,CallABM
USE APHchl_Module
IMPLICIT NONE
CHARACTER(LEN=21), PARAMETER:: ProcName='Combine'
CHARACTER(LEN=10) hour
CHARACTER(LEN=8) today
CHARACTER(LEN=5) curzone
CHARACTER(LEN=3) SFunType(2)
INTEGER NType, NCurves, dtvalues(8)

WRITE(Out_Unit,*)'Called    Combine'
CLOSE(Out_Unit)
OPEN(Unit=Out_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'Output/Combine.txt',Form='FORMATTED')
WRITE(Out_Unit,*)
WRITE(Out_Unit,*)
WRITE(Out_Unit,*)
WRITE(Out_Unit,*)'Entering:', ProcName

CALL Date_And_Time(today, hour, curzone, dtvalues)
CALL Print_Date_And_Time(Today, Hour, CurZone, DtValues, Out_Unit)

IF(CallDVR)THEN
   SFunType(1)='DVR'
ELSEIF(CallFEM)THEN
   SFunType(1)='FEM'
ELSE
   SFunType(1)='ABM'
ENDIF

SFunType(2)='ABM'
IF(SFunType(1)==SFunType(2))THEN
   NType=1
ELSE
   NType=2
ENDIF

! Combine the surface function energy levels into a single file
CALL Combine_SFLevels(SFunType,NType)
! Combine Overlap Matrices into a single file
CALL Combine_Overlap(CallFEM,CallDVR,CallPDAF,CallABM, NCurves)
! Combine Potential Matrices into a single file
CALL Combine_PotMatrix(CallFEM,CallDVR,CallPDAF,CallABM)

WRITE(Out_Unit,*)'Leaving:', ProcName

CALL Date_And_Time(today, hour, curzone, dtvalues)
CALL Print_Date_And_Time(Today, Hour, CurZone, DtValues, Out_Unit)
WRITE(Out_Unit,*)
WRITE(Out_Unit,*)'Completed Combine'
WRITE(Out_Unit,*)

CLOSE(Out_Unit)
OPEN(Unit=Out_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'Output/BasisForScatt.txt',Form='FORMATTED',ACCESS='APPEND')
WRITE(Out_Unit,*)'Completed Combine'
WRITE(Out_Unit,*)

RETURN
ENDSUBROUTINE Combine
