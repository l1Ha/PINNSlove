MODULE DelAng_Module
   USE Numeric_Kinds_Module
   IMPLICIT NONE
   SAVE
   REAL(Kind=WP_Kind) dthe     !  theta angle increment.
   REAL(Kind=WP_Kind) dchi(3)  !  chi angle increment for each channel.
END MODULE DelAng_Module
