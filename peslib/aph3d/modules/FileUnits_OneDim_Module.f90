MODULE FileUnits_OneDim_Module
   USE Numeric_Kinds_Module
   IMPLICIT NONE
   SAVE
   !CHARACTER(LEN=200) InputDir
   !CHARACTER(LEN=200) OutDIr
   INTEGER(KIND=IW_Kind), PARAMETER:: Screen_Unit            = 0
   INTEGER(KIND=IW_Kind), PARAMETER:: Dbug_Unit              = 7
   INTEGER(KIND=IW_Kind), PARAMETER:: Time_Unit              = 8
   INTEGER(KIND=IW_Kind), PARAMETER:: Pot_Unit               =10
   INTEGER(KIND=IW_Kind), PARAMETER:: SYSPot_Unit            =11
   INTEGER(KIND=IW_Kind), PARAMETER:: Crit_Unit              =12
   INTEGER(KIND=IW_Kind), PARAMETER:: Graf_Unit              =13
   INTEGER(KIND=IW_Kind), PARAMETER:: Morse_Unit             =14
   INTEGER(KIND=IW_Kind), PARAMETER:: DAF_Unit               =15
   INTEGER(KIND=IW_Kind), PARAMETER:: DAFBin_Unit            =16
   INTEGER(KIND=IW_Kind), PARAMETER:: Numerov_Unit           =17
   INTEGER(KIND=IW_Kind), PARAMETER:: NumerovBin_Unit        =18
   INTEGER(KIND=IW_Kind), PARAMETER:: Spectra_Unit           =19
   INTEGER(KIND=IW_Kind), PARAMETER:: PsiDAFG_Unit           =20
   INTEGER(KIND=IW_Kind), PARAMETER:: PsiDAFGS_Unit          =21
   INTEGER(KIND=IW_Kind), PARAMETER:: PsiDAFD_Unit           =22
   INTEGER(KIND=IW_Kind), PARAMETER:: PsiDAFDS_Unit          =23
   INTEGER(KIND=IW_Kind), PARAMETER:: PsiNumerovG_Unit       =24
   INTEGER(KIND=IW_Kind), PARAMETER:: PsiNumerovGS_Unit      =25
   INTEGER(KIND=IW_Kind), PARAMETER:: PsiNumerovD_Unit       =26
   INTEGER(KIND=IW_Kind), PARAMETER:: PsiNumerovDS_Unit      =27
   INTEGER(KIND=IW_Kind), PARAMETER:: DAF_DVR_Unit           =28
   INTEGER(KIND=IW_Kind), PARAMETER:: Numerov_DVR_Unit       =29
   INTEGER(KIND=IW_Kind), PARAMETER:: EigPltNumerov_Unit     =30
   INTEGER(KIND=IW_Kind), PARAMETER:: EigPltDAF_Unit         =31
   INTEGER(KIND=IW_Kind), PARAMETER:: Numerov_Prop_Unit      =32
   INTEGER(KIND=IW_Kind), PARAMETER:: Continuum_Unit         =33
   INTEGER(KIND=IW_Kind), PARAMETER:: ContinuumBin_Unit      =34
   INTEGER(KIND=IW_Kind), PARAMETER:: BC_Plot_Unit           =35
   INTEGER(KIND=IW_Kind), PARAMETER:: Eig_Plot_Unit          =36 
   INTEGER(KIND=IW_Kind), PARAMETER:: Dipole_Unit            =37
   INTEGER(KIND=IW_Kind), PARAMETER:: DAF_Oscil_Unit         =38
   INTEGER(KIND=IW_Kind), PARAMETER:: Numerov_Oscil_Unit     =39
   INTEGER(KIND=IW_Kind), PARAMETER:: BC_Bin_Unit            =40
   INTEGER(KIND=IW_Kind), PARAMETER:: Diatomic_Out_Unit      =2966
   INTEGER(KIND=IW_Kind), PARAMETER:: FourierGrid_Unit       =2967
   INTEGER(KIND=IW_Kind), PARAMETER:: FourierGridBin_Unit    =2968
   INTEGER(KIND=IW_Kind), PARAMETER:: EigPltFourierGrid_Unit =2969
   INTEGER(KIND=IW_Kind), PARAMETER:: PsiFourierGridGS_Unit  =2970
   INTEGER(KIND=IW_Kind), PARAMETER:: PsiFourierGridG_Unit   =2971
   INTEGER(KIND=IW_Kind), PARAMETER:: PsiFourierGridDS_Unit  =2972
   INTEGER(KIND=IW_Kind), PARAMETER:: PsiFourierGridD_Unit   =2973
   INTEGER(KIND=IW_Kind), PARAMETER:: FourierGrid_DVR_Unit   =2974
   INTEGER(KIND=IW_Kind), PARAMETER:: FourierGrid_Oscil_Unit =2975
   INTEGER(KIND=IW_Kind), PARAMETER:: YfitPltDaf_Unit        =2976
   INTEGER(KIND=IW_Kind), PARAMETER:: YfitPltNumerov_Unit    =2977
   INTEGER(KIND=IW_Kind), PARAMETER:: YfitPltFourier_Unit    =2978
   INTEGER(KIND=IW_Kind), PARAMETER:: YfitBin_Unit           =2979
   INTEGER(KIND=IW_Kind), PARAMETER:: Yfitter_Unit           =3000
   
END MODULE FileUnits_OneDim_Module