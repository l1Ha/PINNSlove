MODULE Vector_Module
   USE Numeric_Kinds_Module
   SAVE
   INTEGER(KIND=IW_Kind) igr(64)   !
   REAL(Kind=WP_Kind) fpn(64)   !
   NAMELIST/vector/ fpn, igr
END MODULE Vector_Module
