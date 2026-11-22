SUBROUTINE BasisRegions(CallFEM, CallDVR, CallPDAF, CALLABM, CallDelves)
USE Numeric_Kinds_Module
USE Numbers_Module
USE BasisRegions_Module
USE Boundary_Module
USE FileUnits_Module
USE SfnDis_Module
IMPLICIT NONE
!  This routine sets up the propagation regions.
!
!  DVR/FEM/PDAF          ABM_APH         ABM_Delves/Delves             Jacobi
!
!  Distance increases to the right  ----->
!--------------------------------------------------------------------
INTEGER IthSector
LOGICAL CallFEM, CallDVR, CallPDAF, CallABM, CallDelves

WRITE(Out_Unit,*)'Called    BasisRegions'
CLOSE(Out_Unit)
OPEN(Unit=Out_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'Output/BasisRegions.txt',Form='FORMATTED')

WRITE(Out_Unit,*)
WRITE(Out_Unit,'("RhoStart        =",ES23.15)')RhoStart
WRITE(Out_Unit,'("RhoStartABM_APH =",ES23.15)')RhoStartABM_APH
WRITE(Out_Unit,'("RhoStartDelves  =",ES23.15)')RhoStartDelves
WRITE(Out_Unit,'("RStartJacobi    =",ES23.15)')RStartJacobi

WRITE(Out_Unit,*)
WRITE(Out_Unit,'("JacobiProp     =",L2     )')JacobiProp

IF(RhoStartDelves>RStartJacobi)THEN
   WRITE(Out_Unit,'("RhoStartDelves>RStartJacobi",2ES23.15)')RhoStartDelves,RStartJacobi
   WRITE(Msg_Unit,'("RhoStartDelves>RStartJacobi",2ES23.15)')RhoStartDelves,RStartJacobi
   STOP "Stopping in BasisRegions"
ENDIF

IF(RhoStartABM_APH>RhoStartDelves)THEN
   WRITE(Out_Unit,'("RhoStartABM_APH>RhoStartDelves",2ES23.15)')RhoStartABM_APH,RhoStartDelves
   WRITE(Msg_Unit,'("RhoStartABM_APH>RhoStartDelves",2ES23.15)')RhoStartABM_APH,RhoStartDelves
   STOP "Stopping in BasisRegions"
ENDIF

IF(RhoStart>RhoStartABM_APH)THEN
   WRITE(Out_Unit,'("RhoStart>RhoStartABM_APH",2ES23.15)')RhoStart,RhoStartABM_APH
   WRITE(Msg_Unit,'("RhoStart>RhoStartABM_APH",2ES23.15)')RhoStart,RhoStartABM_APH
   STOP "Stopping in BasisRegions"
ENDIF
!--------------------------------------------------------------------
!  DVR/FEM/PDAF Region
!--------------------------------------------------------------------
IF(RhoStart < RhoStartABM_APH)THEN
!NStartSectors=Nsectors !! TMPMOD
   NStartSectors=NINT(NSectors*(RhoStartABM_APH-RhoStart)/(RStartJacobi-RhoStart))
   !tmpmod GregParkerIF(RhoStartABM_APH==RStartJacobi)NStartSectors=NSectors+1
   !IF(NStartSectors==Nsectors)NStartSectors=NStartSectors-1
   NStartSectors=MAX(NStartSectors,2)
   WRITE(Out_Unit,'("DVRBasis       =",L2)')DVRBasis
   WRITE(Out_Unit,'("PDAFBasis      =",L2)')PDAFBasis
   WRITE(Out_Unit,'("FEMBasis       =",L2)')FEMBasis
   WRITE(Out_Unit,*)
   IF(DVRBasis)THEN
      WRITE(Out_Unit,'("DVR_Basis starts at",ES23.15, " and ends at",ES23.15)')RhoStart,RhoStartABM_APH
      WRITE(Out_Unit,'("The DVR_to_ABM Transformation is generated at" ,ES23.15)')RhoStartABM_APH
   ELSEIF(PDAFBasis)THEN
      WRITE(Out_Unit,'("PDAF_Basis starts at",ES23.15," and ends at",ES23.15)')RhoStart,RhoStartABM_APH
      WRITE(Out_Unit,'("The PDAF_to_ABM Transformation is generated at",ES23.15)')RhoStartABM_APH
   ELSEIF(FEMBasis)THEN
      WRITE(Out_Unit,'("FEM_Basis starts at",ES23.15, " and ends at",ES23.15)')RhoStart,RhoStartABM_APH
      WRITE(Out_Unit,'("The FEM_to_ABM Transformation is generated at" ,ES23.15)')RhoStartABM_APH
   ELSE
      NStartSectors=0
      WRITE(Out_Unit,'("Error: DVRBasis, FEMBasis and PDAFBasis are all False")')
      WRITE(Msg_Unit,'("Error: DVRBasis, FEMBasis and PDAFBasis are all False")')
      STOP 'Stopping in BasisRegions'
   ENDIF
   WRITE(Out_Unit,'("NStartSectors  =",I5)')NStartSectors
   CALL Boundaries(RhoStart, RhoStartABM_APH, DeltaRho1, DeltaRho2, NStartSectors)
   DO IthSector=1,NStartSectors
      IF(FEMBasis) Basis_Type(IthSector)='FEM'
      IF(PDAFBasis)Basis_Type(IthSector)='PDAF'
      IF(DVRBasis) Basis_Type(IthSector)='DVR'
   ENDDO
   IF(FEMBasis)THEN
      IthStart_FEM=1
      IthEnd_FEM=NStartSectors
   ENDIF
   IF(DVRBasis)THEN
      IthStart_DVR=1
      IthEnd_DVR=NStartSectors
   ENDIF
   Sector_Boundaries(NStartSectors+1)=MIN(RhoStartABM_APH,RhoStartDelves)
   Basis_Dist(NStartSectors)=(Sector_Boundaries(NStartSectors)+Sector_Boundaries(NStartSectors+1))/2.d0
ELSE
   WRITE(Out_Unit,'("NO DVR/FEM/PDAF Region!")')
ENDIF
WRITE(Out_Unit,*)

!--------------------------------------------------------------------
!  APH_ABM Region
!--------------------------------------------------------------------
IF(RhoStartABM_APH <= RhoStartDelves)THEN
   IthStart_ABM=MAX(IthEnd_FEM,IthEnd_DVR)+1
   NAPH_ABMSectors=NINT(NSectors*(RhoStartDelves-RhoStartABM_APH)/(RStartJacobi-RhoStart))
   IF(.Not.(RhoStartDelves< RStartJacobi))NAPH_ABMSectors=NSectors-NStartSectors
   NAPH_ABMSectors=MAX(NAPH_ABMSectors,2)
   WRITE(Out_Unit,'("APH_ABM_Basis starts at",ES23.15," and ends at",ES23.15)')RhoStartABM_APH,RhoStartDelves
   WRITE(Out_Unit,'("The ABM_to_Delves transformation is generated at",ES23.15)')RhoStartDelves
   WRITE(Out_Unit,'("NAPH_ABMSectors    =",I5)')NAPH_ABMSectors
   CALL Boundaries(RhoStartABM_APH, RhoStartDelves, DeltaRho1, DeltaRho2, NAPH_ABMSectors)
   DO IthSector=NStartSectors+1,NStartSectors+NAPH_ABMSectors
      Basis_Type(IthSector)='ABM_APH'
   ENDDO
   IthStart_ABM=NStartSectors+1
   IthEnd_ABM=NStartSectors+NAPH_ABMSectors
ELSE
   NAPH_ABMSectors=0
   WRITE(Out_Unit,'("NO APH Region!")')
ENDIF
WRITE(Out_Unit,*)
IF(IthEnd_FEM>=2.or.IthEnd_DVR>=2)THEN
   IF(IthStart_ABM>=2)THEN
      Basis_Dist(IthStart_ABM)=(Sector_Boundaries(IthStart_ABM)+Sector_Boundaries(IthStart_ABM+1))/2.d0
   ELSEIF(IthEnd_DVR==1)THEN
      Basis_Dist(IthStart_DVR)=(Sector_Boundaries(IthStart_DVR)+Sector_Boundaries(IthStart_DVR+1))/2.d0
   ELSEIF(IthEnd_FEM==1)THEN
      Basis_Dist(IthStart_FEM)=(Sector_Boundaries(IthStart_FEM)+Sector_Boundaries(IthStart_FEM+1))/2.d0
   ELSEIF(IthEnd_PDAF==1)THEN
      Basis_Dist(IthStart_PDAF)=(Sector_Boundaries(IthStart_PDAF)+Sector_Boundaries(IthStart_PDAF+1))/2.d0
   ENDIF
ENDIF

!--------------------------------------------------------------------
!  Delves Region
!--------------------------------------------------------------------
IF(RhoStartDelves < RStartJacobi)THEN
   NDelvesSectors=NSectors-NStartSectors-NAPH_ABMSectors+1
   WRITE(Out_Unit,'("Delves_Basis starts at",ES23.15," and ends at",ES23.15)')RhoStartDelves,RStartJacobi
   WRITE(Out_Unit,'("The Delves_to_Jacobi Transformation is generated at",ES23.15)')RStartJacobi
   WRITE(Out_Unit,'("NDelvesSectors =",I5)')NDelvesSectors
   CALL Boundaries(RhoStartDelves, RStartJacobi, DeltaRho1, DeltaRho2, NDelvesSectors)
   DO IthSector=NStartSectors+NAPH_ABMSectors+1,NStartSectors+NAPH_ABMSectors+NDelvesSectors
      Basis_Type(IthSector)='Delves'
   ENDDO
ELSE
   NDelvesSectors=0
   WRITE(Out_Unit,'("NO Delves Region!")')
ENDIF
WRITE(Out_Unit,*)

!--------------------------------------------------------------------
!  Jacobi Region
!--------------------------------------------------------------------
WRITE(Out_Unit,*)
WRITE(Out_Unit,'("JacobiProp     =",L2     )')JacobiProp
IF(JacobiProp)THEN
   WRITE(Out_Unit,'("Jacobi_Basis starts at",ES23.15," and ends at Infinity")')RStartJacobi
ELSE
   WRITE(Out_Unit,'("NO Jacobi Region!")')
ENDIF
WRITE(Out_Unit,*)

CallFEM    =.False.
CallDVR    =.False.
CallPDAF   =.False.
CallABM    =.False.
CallDelves =.False.
IF(RhoStart < RStartJacobi)THEN
   WRITE(Out_Unit,'("Basis Functions will be calculated at NSectors =",I5)')NSectors
   Basis_Dist(NSectors)=RStartJacobi
   DO IthSector=1,NSectors
      WRITE(Out_Unit,'("IthSector=",I5," BasisDist=",ES22.15," Basis_Type=",A10)') &
            IthSector,Basis_Dist(IthSector),    &
            Basis_Type(IthSector)
      IF(Basis_Type(IthSector)=='FEM'    )CallFEM   =.True.
      IF(Basis_Type(IthSector)=='DVR'    )CallDVR   =.True.
      IF(Basis_Type(IthSector)=='PDAF'   )CallPDAF  =.True.
      IF(Basis_Type(IthSector)=='ABM_APH')CallABM   =.True.
      IF(Basis_Type(IthSector)=='Delves' )CallDelves=.True.
   ENDDO
   WRITE(Out_Unit,*)
   WRITE(Out_Unit,'("The Sector Boundaries are:")')
   IF(JacobiProp)Sector_Boundaries(NSectors+1)=RStartJacobi
   IF(JacobiProp)Sector_Boundaries(NSectors+1)=Ten**(RANGE(One)-10)
   DO IthSector=1,NSectors
      WRITE(Out_Unit,'("IthSector=",I5," Left Sector Boundary=",ES22.15," Right Sector Boundary=",ES22.15)')   &
                        IthSector,Sector_Boundaries(IthSector),Sector_Boundaries(IthSector+1)
   ENDDO
   WRITE(Out_Unit,*)
   WRITE(Out_Unit,'("Adiabatic curves start at",ES23.15," and end at",ES23.15)')RhoStart,RStartJacobi
   WRITE(Out_Unit,'("Diabatic curves start at ",ES23.15," and end at",ES23.15)')RhoStart,RStartJacobi
   WRITE(Out_unit,'("A 1-Dimensional Analysis of the Curves will be performed")')
ENDIF

CLOSE(Out_Unit)
OPEN(Unit=Out_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'Output/BasisForScatt.txt',Form='FORMATTED',ACCESS='APPEND')
WRITE(Out_Unit,*)'Completed BasisRegions'
WRITE(Out_Unit,*)

RETURN
ENDSUBROUTINE BasisRegions
