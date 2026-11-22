MODULE DBS_Module
   USE Numeric_Kinds_Module
   SAVE
   REAL(Kind=WP_Kind) ksqmax    !
   INTEGER(KIND=IW_Kind) nsteps ! Number of steps per asymptotic wavelength used in the Log-derivative method.
   LOGICAL frst_WRITE           !
END MODULE DBS_Module
