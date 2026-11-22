MODULE Cook_Module
   USE Numeric_Kinds_Module
   SAVE
   REAL(Kind=WP_Kind) ycord   !
   REAL(Kind=WP_Kind) zcord   !
   REAL(Kind=WP_Kind) enint   !
   INTEGER(KIND=IW_Kind) icook   !
   NAMELIST/cook/ enint, ycord, zcord, icook
END MODULE Cook_Module
