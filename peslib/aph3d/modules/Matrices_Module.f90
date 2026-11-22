MODULE Matrices_Module
   USE Numeric_Kinds_Module
   USE Max_Module
   IMPLICIT NONE
   SAVE
   REAL(Kind=WP_Kind) rd(maxang,maxang)   !
   REAL(Kind=WP_Kind) rdt(maxang,maxang)  !
   REAL(Kind=WP_Kind) scon(maxang,maxang) !
   REAL(Kind=WP_Kind) aa(maxang,maxang)   !
   REAL(Kind=WP_Kind) subd(maxang)        !
   REAL(Kind=WP_Kind) sd(maxangv)         !
   REAL(Kind=WP_Kind) en(maxang)          !
END MODULE Matrices_Module
