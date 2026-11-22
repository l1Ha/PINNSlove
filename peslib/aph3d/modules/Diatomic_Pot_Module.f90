MODULE DiatomicPot_MODULE
   USE Numeric_Kinds_Module
   USE Narran_Module
   IMPLICIT NONE
   SAVE
   CHARACTER(LEN=25) DiatomicPot   ! A unique name for the potential energy.
   LOGICAL FixedMaximum
   INTEGER(KIND=IW_KIND) NStatesDAF(Narran), NstatesNumerov(Narran), nstatesfouriergrid(Narran)
   REAL(KIND=WP_KIND) FixedMax
END MODULE DiatomicPot_MODULE