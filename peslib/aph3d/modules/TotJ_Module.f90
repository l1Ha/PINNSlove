MODULE TotJ_Module
   USE Narran_Module
   USE Numeric_Kinds_Module
   USE CSBasis_Module
   IMPLICIT NONE
   SAVE
   LOGICAL:: Symmetry=.TRUE.  ! True if symmetry decoupling is being used.  This can be the case only if channels B and C are identical atoms
   LOGICAL:: jEven=.TRUE.     ! Ignored if symmetry is false.
                                      ! If symmetry is true and jEven is true then even rotational states are used in channel A.
                                      ! If symmetry is true and jEven is false then odd rotational states are used in channel A.
   LOGICAL:: Nsymc(Narran)
   INTEGER(KIND=IW_Kind):: JTot=0           ! Total angular momentum, J.
   INTEGER(KIND=IW_Kind):: Parity=0         ! Parity.
   INTEGER(KIND=IW_Kind):: Mega=0           ! Total angular momentum projection quantum number.
   INTEGER(KIND=IW_Kind):: JRef=0           ! Reference J used to get surface functions independent of the actual J.
   INTEGER(KIND=IW_Kind):: MegaCoup=0
   INTEGER(KIND=IW_Kind):: MegaMin=0        ! Minimum value of the projection quantum number.
   INTEGER(KIND=IW_Kind):: MegaMax=0        ! Maximum value of the projection quantum number.
   NAMELIST/Momentum/ JTot, Parity, jEven, MegaMax, NSymc, JRef, MegaCoup
   NAMELIST/totjb/ jtot, parity, mega, megamin, megamax, symmetry, jeven, jref, megacoup
END MODULE TotJ_Module
