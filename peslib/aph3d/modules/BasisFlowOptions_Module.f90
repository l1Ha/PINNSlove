MODULE BasisFlowOptions_Module
   USE Numeric_Kinds_Module
   USE BasisRegions_Module
   USE SfnDis_Module
   LOGICAL :: CallPlotVee    = .False.
   LOGICAL :: CallBound3D    = .False.
   LOGICAL :: CalcSFunc      = .True.
   LOGICAL :: CallFEM        = .False.
   LOGICAL :: CallDVR        = .False.
   LOGICAL :: CallPDAF       = .False.
   LOGICAL :: CallABM        = .True.
   LOGICAL :: CallDelves     = .False.
   LOGICAL :: CallPropagate  = .True.
   LOGICAL :: CallAsymptotic = .True.
   LOGICAL :: CallDiatomic   = .True.
   LOGICAL :: There
   INTEGER :: ISector
   INTEGER :: dtvalues(8)
   ! PORT NOTE (2026-10-01): CallBound3D also listed in the SectorDef group;
   ! the shipped input decks (e.g. APH3D_H3_DMBE/DMBE_ABM.data) put it there.
   NAMELIST/SectorDef/   RhoStart, RhoStartABM_APH, RhoStartDelves,RStartJacobi, NSectors, JacobiProp, DVRBasis, FEMBasis, PDAFBasis, Kase_CalcParm, DeltaRho1, DeltaRho2, CallBound3D
   NAMELIST/FlowOptions/ CallPlotVee, CallBound3D, CalcSFunc, CallFEM, CallDVR, CallPDAF, CallABM, CallDelves, CallPropagate, CallAsymptotic, CallDiatomic
   !NAMELIST/sfundist/ rhomin, rhomax, deltarho1, deltarho2, nrho See SfnDis_Module
   !rhomin=RhoStart
   !rhomax=RStartJacobi
   !DeltaRho1 (Namelist read)
   !DeltaRho2 (Namelist read)
   !nrho=NSectors
END MODULE BasisFlowOptions_Module