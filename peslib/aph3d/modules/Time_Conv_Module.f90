MODULE Time_Conv_Module
   USE Numeric_Kinds_Module
   IMPLICIT NONE
   SAVE
   REAL(Kind=WP_Kind) time_conv   !
   REAL(Kind=WP_Kind) hbar        !
   REAL(Kind=WP_Kind) hartree     !
   NAMELIST/timefac/ time_conv  !
END MODULE Time_Conv_Module
