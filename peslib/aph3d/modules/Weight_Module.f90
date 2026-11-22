MODULE Weight_Module
   USE Narran_Module
   USE Numeric_Kinds_Module
   SAVE
   REAL(Kind=WP_Kind) betaw(narran)  !
   REAL(Kind=WP_Kind) seq(narran)    !
   REAL(Kind=WP_Kind) chncof(narran) !
   NAMELIST/weight/ betaw, seq, chncof
END MODULE Weight_Module
