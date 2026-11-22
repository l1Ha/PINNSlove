MODULE MassFactor2_Module
   USE Numeric_Kinds_Module
   USE Parms_Module
   IMPLICIT NONE
   SAVE
   INTEGER(KIND=IW_Kind) jacute(narran)       !
   REAL(Kind=WP_Kind) d(narran)            !
   REAL(Kind=WP_Kind) acute(narran)        !
   REAL(Kind=WP_Kind) r(narran)            !
   REAL(Kind=WP_Kind) chij(narran, narran) !
   REAL(Kind=WP_Kind) s(narran)            !
   REAL(Kind=WP_Kind) chi3(narran)         !
END MODULE MassFactor2_Module