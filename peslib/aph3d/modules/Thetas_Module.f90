MODULE Thetas_Module
   USE Narran_Module
   USE Numeric_Kinds_Module
   SAVE
   REAL(Kind=WP_Kind) c(narran)  !
   REAL(Kind=WP_Kind) b(narran)  !
   NAMELIST/thetas/ c,b
END MODULE Thetas_Module
