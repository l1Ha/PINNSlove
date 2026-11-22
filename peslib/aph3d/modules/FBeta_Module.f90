MODULE FBeta_Module
   USE Numeric_Kinds_Module
   USE Narran_Module
   SAVE 
   ! Array containing the beta values used at the first rho in
   ! factoring out the exponential behaivor of the surface
   ! functions.
   REAL(Kind=WP_Kind) :: betafrst(narran) = [10.0D0,10.0D0,10.0D0]
   NAMELIST/fbeta/ betafrst
END MODULE FBeta_Module
