MODULE PotZero_Module
   USE Numeric_Kinds_Module
   SAVE
   REAL(Kind=WP_Kind):: vzero=0.d0  ! Amount to shift the potential
   NAMELIST/potzero/ vzero
END MODULE PotZero_Module
