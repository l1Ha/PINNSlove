MODULE FileUnits_Asymptotic_Module
   USE Numeric_Kinds_Module
   USE FileUnits_Common_Module
   IMPLICIT NONE
   SAVE
   INTEGER(KIND=IW_Kind), PARAMETER:: DerUnit          =  87        !
   INTEGER(KIND=IW_Kind), PARAMETER:: EigUnit          =  88        !
   INTEGER(KIND=IW_Kind), PARAMETER:: Prob_Unit        = 103        !
   INTEGER(KIND=IW_Kind), PARAMETER:: Probf_Unit       = 104        !
   INTEGER(KIND=IW_Kind), PARAMETER:: KSym_Unit        = 301
   INTEGER(KIND=IW_Kind), PARAMETER:: KESym_Unit       = 302
   INTEGER(KIND=IW_Kind), PARAMETER:: SReal_Unit       = 303
   INTEGER(KIND=IW_Kind), PARAMETER:: SImag_Unit       = 304
   INTEGER(KIND=IW_Kind), PARAMETER:: SEReal_Unit      = 305
   INTEGER(KIND=IW_Kind), PARAMETER:: SEImag_Unit      = 306
   INTEGER(KIND=IW_Kind), PARAMETER:: Probil_Unit      = 307
   INTEGER(KIND=IW_Kind), PARAMETER:: ProbilE_Unit     = 308
   INTEGER(KIND=IW_Kind), PARAMETER:: Phase_Unit       = 309
   INTEGER(KIND=IW_Kind), PARAMETER:: PhaseE_Unit      = 310
   INTEGER(KIND=IW_Kind), PARAMETER:: PhaseTran_Unit   = 311
   INTEGER(KIND=IW_Kind), PARAMETER:: PhaseTranE_Unit  = 312
   INTEGER(KIND=IW_Kind), PARAMETER:: QReal_Unit       = 313
   INTEGER(KIND=IW_Kind), PARAMETER:: QImag_Unit       = 314
   INTEGER(KIND=IW_Kind), PARAMETER:: TDelayReal_Unit  = 315
   INTEGER(KIND=IW_Kind), PARAMETER:: TDelayImag_Unit  = 316
   INTEGER(KIND=IW_Kind), PARAMETER:: EigLifePltUnit   = 318
   INTEGER(KIND=IW_Kind), PARAMETER:: PhasePltUnit     = 319
   INTEGER(KIND=IW_Kind), PARAMETER:: PhaseEPltUnit    = 320
   INTEGER(KIND=IW_Kind), PARAMETER:: Plt_Unit         = 401
   INTEGER(KIND=IW_Kind), PARAMETER:: AP_Unit          = 402
   INTEGER(KIND=IW_Kind), PARAMETER:: In_UnitA         = 403
   INTEGER(KIND=IW_Kind), PARAMETER:: In_UnitB         = 404
END MODULE FileUnits_Asymptotic_Module