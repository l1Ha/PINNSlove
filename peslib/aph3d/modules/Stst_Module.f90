MODULE  Stst_Module
   USE Numeric_Kinds_Module
   USE Narran_Module
   USE Quad_Module
   USE FBeta_Module
   IMPLICIT NONE
   SAVE
   REAL(Kind=WP_Kind) :: Seq(NArran)    = [1.13532D0,1.13532D0,1.13532D0] !
   REAL(Kind=WP_Kind) :: Betaw(NArran)  = [21.51D0,21.51000D0,21.51D0]    !
   REAL(Kind=WP_Kind) :: chncof(NArran) = [1.0D0,1.0D0,1.0D0]             !
   NAMELIST/stst/ nquad, seq, betaw, chncof, betafrst
END MODULE Stst_Module
