MODULE Integopt_Module
   USE Numeric_Kinds_Module
   IMPLICIT NONE
   SAVE
   CHARACTER(LEN=13) int_opt   !
   REAL(Kind=WP_Kind) rhodvr   !
   REAL(Kind=WP_Kind) rhoabm   !
   REAL(Kind=WP_Kind) rhofem   !
   NAMELIST/integopt/ rhodvr, rhoabm, rhofem, int_opt
END MODULE Integopt_Module
