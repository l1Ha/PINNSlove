INTEGER FUNCTION CalcParm(NRho, Start, DeltaRho1, DeltaRho2, Finish, Kase_CalcParm)
USE Numeric_Kinds_Module
USE Numbers_Module
USE FileUnits_Module
IMPLICIT NONE
SAVE
CHARACTER(LEN=21), PARAMETER:: ProcName='CalcParm'
REAL(Kind=WP_Kind) Start, DeltaRho1, DeltaRho2, Finish
INTEGER NRho, Kase_CalcParm

WRITE(Out_Unit,*)
WRITE(Out_Unit,*)'   CalcParm: Start'

CalcParm = 0
IF(NRho==1)THEN
   IF(Kase_CalcParm==0)Finish=Start
   IF(Kase_CalcParm==1)DeltaRho1 = Zero
   IF(Kase_CalcParm==2)DeltaRho2 = Zero
   CalcParm=0
ELSE
   CalcParm=0
   IF(Kase_CalcParm==0)THEN
      Finish = (Start + (NRho - 1)*DeltaRho1)* (1 + DeltaRho2)**(NRho - 1)
   ELSEIF(Kase_CalcParm==1)THEN
      DeltaRho1 = (Finish/(1 + DeltaRho2)**(NRho - 1) - Start) / (NRho - 1)
   ELSEIF(Kase_CalcParm==2)THEN
      DeltaRho2 = (Finish/(Start + (NRho - 1)*DeltaRho1)) ** (1.d0/(NRho-1)) - One
   ELSE
      WRITE(Msg_Unit,*) 'Error Kase_CalcParm incorrect'
      CalcParm = -1
   ENDIF
ENDIF

WRITE(Out_Unit,'("   NRho=",I5," Start=",ES23.15," Finish=",ES23.15," Kase_CalcParm=",I2)')NRho, Start, Finish, Kase_CalcParm
WRITE(Out_Unit,'("   DeltaRho1=",ES23.15," DeltaRho2=",ES23.15)')DeltaRho1, DeltaRho2
WRITE(Out_Unit,*)'   CalcParm: END'
WRITE(Out_Unit,*)

IF(DeltaRho1<Zero.or.DeltaRho2<Zero)THEN
   WRITE(*,*)"Error for input values in CalcPark"
   WRITE(*,*)"DeltaRho1 and DeltaRho2 must be great than Zero", DeltaRho1, DeltaRho2
   WRITE(*,*)"Error for input values in CalcPark"
   WRITE(*,*)"DeltaRho1 and DeltaRho2 must be great than Zero", DeltaRho1, DeltaRho2
   STOP "Error in CalcParm"
ENDIF

RETURN
END FUNCTION CalcParm
