SUBROUTINE BasisForScattering
USE Numeric_Kinds_Module
USE Regins_Module   ! PORT (2026-10-02): endaph/enddelve
USE BasisRegions_Module
USE Boundary_Module
USE FileUnits_Module
USE PES_MODULE
USE InputFile_Module
USE Masses_Module
USE NumNuj_Module, ONLY: nnuj
USE GaussB_Module, ONLY: NHermt
USE BasisFlowOptions_Module
USE Region_Module
USE Oops_Module
USE APHChl_Module, ONLY : neigmin_sum, naph
IMPLICIT NONE
CHARACTER(LEN=10) hour
CHARACTER(LEN=8) today
CHARACTER(LEN=5) curzone
INTEGER NCurves, NopenMax

OPEN(Unit=DBug_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'Output/debuging.txt',Form='FORMATTED')
OPEN(Unit=Out_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'Output/BasisForScatt.txt',Form='FORMATTED')
CALL Date_And_Time(today, hour, curzone, dtvalues)                  ! Write Time and Data
WRITE(Out_Unit,'("At the start of AAA_Basis_Main")')
CALL Print_Date_And_Time(Today, Hour, CurZone, DtValues, Out_Unit)

WRITE(Out_Unit,*)'Read in All of the Input Data'
CALL ReadAllData
WRITE(*,*) 'DBG BasisForScattering after ReadAllData: finish=', finish, ' endaph=', endaph
WRITE(Out_Unit,*) 'DBG BasisForScattering after ReadAllData: finish=', finish, ' endaph=', endaph

WRITE(Out_Unit,*)'Determine Basis Set Distances'
CALL BasisRegions(CallFEM, CallDVR, CallPDAF, CallABM, CallDelves)
  
IF(CallPlotVee)THEN 
   WRITE(*,*)"Starting CallPlotVee Block"
   WRITE(Out_Unit,*)"Starting CallPlotVee Block"
   
   WRITE(Out_Unit,*)'Plot Collinear Potential at each sector'
   WRITE(Msg_Unit,*)'Plot Collinear Potential at each sector'
   CALL CollinearFixedRho(Basis_Dist, NSectors)
   
   WRITE(Out_Unit,*)'Plot Collinear Potential'
   WRITE(Msg_Unit,*)'Plot Collinear Potential'
   CALL Graf_Collinear(Basis_Dist(NSectors))
   
   WRITE(Out_Unit,*)'Plot Planar Slices of the Potential '
   WRITE(Msg_Unit,*)'Plot Planar Slices of the Potential '
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
   WRITE(*,*)"End of CallPlotVee Block"
   WRITE(Out_Unit,*)"End of CallPlotVee Block"
ENDIF

IF(CallDiatomic)THEN
   WRITE(*,*)"Starting CallDiatomic Block"
   WRITE(Out_Unit,*)"Calculate Bound States of the Diatomic Fragments"
   CALL Diatomic
   CALL BoundDiatomic
   iregion='delves'
   ioops=.true.
   CALL Delves_Basis
   iregion='jacobi'
   ioops=.false.
   CALL Jacobi_Basis
   CALL CompareMethods
ENDIF
!CALL Delves_Basis_Read
!CALL Jacobi_Basis_Read
!stop "tmp stop in basisforscattering"

IF(CallBound3D)THEN
   WRITE(*,*)"Starting CallBound3D Block"
   WRITE(Out_Unit,*)"Starting CallBound3D Block"
   WRITE(Out_Unit,*)"CallBound3d=",CallBound3D
   Call Bound3D
   WRITE(*,*)"End of CallBound3D Block"
   WRITE(Out_Unit,*)"End of CallBound3D Block"
ENDIF

IF(CalcSFunc)THEN
   WRITE(*,*)"Starting CalcSfunc Block"
   WRITE(Out_Unit,*)"Starting CalcSfunc Block"
   IF(CallFEM)THEN
      WRITE(Out_Unit,*)"FEM Surface Functions"
      CALL FEM
   ELSEIF(CallDVR)THEN
      WRITE(Out_Unit,*)"DVR Surface Functions"
      CALL DVR
      IF(.NOT.callabm)CALL APH_to_Delves
   ELSEIF(CallPDAF)THEN
      WRITE(Out_Unit,*)"PDAF Surface Functions"
      !!!   CALL PDAF   !Not working yet
   ENDIF

   IF(CallABM)THEN
      WRITE(Out_Unit,*)"Calculate ABM_APH Surface Functions"
      CALL ABM_APH
   ENDIF
   
   IF(CallFEM)THEN 
      WRITE(Out_Unit,*)"Calculate FEM to ABM Overlap"
      Call FEM_to_ABM !Not working yet
   ENDIF
      
   IF(CallDVR)THEN 
      WRITE(Out_Unit,*)"Calculate DVR to ABM Overlap"
      CALL DVR_to_ABM
   ENDIF
      
   !!!   IF(CallPDAF)CALL PDAF_to_ABM !Not Working yet
   WRITE(Out_Unit,*)"Calculate APH to Delves Transformation"
   CALL APH_to_Delves
   WRITE(Out_Unit,*)"Calculate Delves Surface Functions"
   CALL Delves_Basis

   IF(CallDelves)THEN
      WRITE(Out_Unit,*)"ABM_Delves Surface Functions"
      !!!   CALL ABM_Delves !Not Working yet
      !CALL DelvProp
   ENDIF

   WRITE(Out_Unit,*)"Combine SFLevels, Ovrlp, PMatrx"
   CALL Combine
   ENDIF

   IF (.NOT. CalcSFunc) THEN
      WRITE(Out_Unit,*)"Recomputing APH to Delves Transformation with NBasisFac=3"
      CALL APH_to_Delves
      CALL Delves_Basis
   ENDIF

   IF(CalcSFunc)THEN
   IF(CallFEM.or.CallDVR.or.CallPDAF.or.CallABM.or.CallDelves)THEN
      WRITE(Out_Unit,*)"Use sflevels to generate plotting files"
      CALL SFLevels('GraphicsOut/Sfelevl_All.rbw')
      WRITE(Out_Unit,*)"Construct Adiabatic Curves"
   ENDIF
   
   WRITE(*,*)"End of CalcSfunc Block"
   WRITE(Out_Unit,*)"End of CalcSfunc Block"
ENDIF

WRITE(Out_Unit,*)"Calculate Asymptotic Jacobi Functions"
CALL Jacobi_Basis

IF(CallPropagate)THEN
   WRITE(*,*)"Starting Propagate Block"
   WRITE(Out_Unit,*)"Starting Propagate Block"
   WRITE(Out_Unit,*)"Generating Propagation Matrices"
   CALL main_matgen

   WRITE(Out_Unit,*)"Propagating Scattering Solution"
   CALL APHlogder

   WRITE(Out_Unit,*)"DO URU transformation to get the R-Matrix in Delves Coordinates"
   CALL UMatph

   ! PORT (2026-10-10): Long-range van der Waals propagation hook
   BLOCK
     CHARACTER(LEN=32) :: env_rf, env_rm
     REAL(Kind=WP_Kind) :: rm_val, rf_val
     INTEGER :: stat_rf, stat_rm
     rm_val = 8.1d0
     rf_val = 0.0d0
     CALL get_environment_variable('APH3D_RF', env_rf, status=stat_rf)
     IF (stat_rf == 0) READ(env_rf, *) rf_val
     CALL get_environment_variable('APH3D_RM', env_rm, status=stat_rm)
     IF (stat_rm == 0) READ(env_rm, *) rm_val
     IF (rf_val > rm_val) THEN
        CALL Jacobi_LongRange_Prop(rm_val, rf_val)
     END IF
   END BLOCK

   WRITE(Out_Unit,*)"Generate K_Matrix"
   CALL Kmatrx_Main(NOpenMax)
   WRITE(*,*)"End of Propagate Block"
   WRITE(Out_Unit,*)"End of Propagate Block"
ELSE
   OPEN(Unit=NOpenMAX_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/NOpenMax.bin', form='unformatted', status='unknown')
   READ(NOpenMAX_Unit)NOpenMax
   !nopenmax=naph  !tmpmod gregparker
   naph=nopenmax
   CLOSE(NOpenMAX_Unit)
ENDIF

IF(CallAsymptotic)THEN
   nopenmax=SUM(nnuj/nhermt)
   WRITE(Out_Unit,*)"Asymptotic Analysis to get S_Matrix: NOpenMax=",NOpenMax
   WRITE(Msg_Unit,*)
   WRITE(Msg_Unit,*)"Asymptotic Analysis to get S_Matrix: NOpenMax=",NOpenMax
   CALL Asymptotic_Analysis(NOpenMax)
   WRITE(Out_Unit,*)"Construct Adiabatic and Diabatic Curves"
   CALL Adiabatic(NSectors)
   CALL AnalysisOfCurves
ENDIF

WRITE(Out_unit,*)"At the end of BasisForScattering"
CALL Date_And_Time(today, hour, curzone, dtvalues)
CALL Print_Date_And_Time(Today, Hour, CurZone, DtValues, Out_Unit)
WRITE(Msg_Unit,*)
WRITE(Msg_Unit,*)"End of Execution"
WRITE(Msg_Unit,*)"Completed BasisForScattering"
WRITE(Out_Unit,*)"End of Execution"
WRITE(Out_Unit,*)"Completed BasisForScattering"

RETURN
ENDSUBROUTINE BasisForScattering
