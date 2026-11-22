MODULE NameList_Module
   USE  APHchl_Module
   USE  Approx_Module
   USE  Arrch_Module
   USE  BasisFlowOptions_Module
   USE  CHerm1_Module
   USE  CHerm2_Module
   USE  Chltot_Module
   USE  Contour_Module
   USE  Converge_Module
   USE  Cook_Module
   USE  CSBasis_Module
   USE  Dip_Module
   USE  Dipole_Module
   USE  DoEnlvls_Module
   USE  DVRNL_Module
   USE  Elemnt_Module
   USE  Energy_Module
   USE  Engpro_Module
   USE  Factrl_Module
   USE  FBeta_Module
   USE  FBR_Module
   USE  GaussB_Module
   USE  InputFile_Module
   USE  Integopt_Module
   USE  Integrat_Module
   USE  Iopt_Module
   USE  LogDer_Module
   USE  Masses_Module
   USE  Melm_Module
   USE  MemAdd_Module
   USE  Mshdat_Module
   USE  New_Module
   USE  NumNuj_Module
   USE  Oops_Module
   USE  Opts_Module
   USE  PBASIS_MODULE
   USE  PES_Module
   USE  PlotNL_Module
   USE  PotZero_Module
   USE  Pow_Module
   USE  Propag_Module
   USE  QCase_Module
   USE  QDiv_Module
   USE  Qlam_Module
   USE  Qnumber_Module
   USE  Qnumsp_Module
   USE  Quad_Module
   USE  Quantb_Module
   USE  Racah_Module
   USE  Regins_Module
   USE  Region_Module
   USE  RhoSur_Module
   USE  RhoValue_Module
   USE  RotMax_Module
   USE  SavNam_Module
   USE  Scalag_Module
   USE  SfnDis_Module
   USE  SFnTyp_Module
   USE  Spectro_Module
   USE  Storage_Module
   USE  Stst_Module
   ! System_Module re-exports DVRMasses' zeta/eta; import only its own names
   ! and the rest of DVRMasses so ReadAllData sees a single eta/zeta.
   USE  System_Module, ONLY: pointg, debug, medium, matelem, h3sys, logspace, &
        peigv, fbrtr, nthesec, nchisec, ltheta, nsfunc, nrow, ncol,           &
        ngood1, ngood2, Num_DVR_Sectors, initsec, nrhomx, nfbr, rendsec,      &
        potmsec, pevmsec, rhofbr1, rhofbr2, fmpotma, deltar, rho1, system,    &
        scale, sine, cosn, tang, rot, rmatch, rstar, dstar, req, redmas,      &
        twomu, qchan, dmatch, rmass, xstar, mass
   USE  Testing_Module
   USE  Thetas_Module
   USE  Time_Conv_Module
   USE  Title1_Module
   USE  Tlogd_Module
   USE  ToDim_Module
   USE  TotalEng_Module
   USE  TotJ_Module
   USE  Var_Module
   USE  Vector_Module
   USE  VFunc_Module
   USE  VIVS_Module
   USE  VWaste_Module
   ! PORT NOTE (2026-10-01): pure re-export hub consumed only by ReadAllData.
   ! Restored public re-export; conflicting duplicate names (eta, zeta, medium,
   ! debug) are only a problem at the point of reference, and none of them are
   ! referenced by ReadAllData.
   SAVE
ENDMODULE NameList_Module
