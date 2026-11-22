MODULE ToDim_Module
   USE Numeric_Kinds_Module
   SAVE
   REAL(Kind=WP_Kind) bet   !
   REAL(Kind=WP_Kind) thic  !
   REAL(Kind=WP_Kind) de    !
   INTEGER(KIND=IW_Kind) iel   !
   INTEGER(KIND=IW_Kind) nnd5  !
   NAMELIST/todim/ bet, thic, de, iel, nnd5
END MODULE ToDim_Module
