MODULE Qnumber_Module
   USE Narran_Module
   USE NState_Module
   USE Numeric_Kinds_Module
   IMPLICIT NONE
   SAVE
   INTEGER(KIND=IW_Kind) :: mvib(nstate)    !  vibrational quantum numbers
   INTEGER(KIND=IW_Kind) :: jthrot(nstate)  ! rotor quantum numbers
   INTEGER(KIND=IW_Kind) :: lorb(nstate)    ! orbital quantum numbers
   INTEGER(KIND=IW_Kind) :: megaz(nstate)   ! projection quantum numbers
   INTEGER(KIND=IW_Kind) :: nuj(nstate)     ! state indices
   INTEGER(KIND=IW_Kind) :: mtype(nstate)   ! Symmetry type for Omegas
   REAL(Kind=WP_Kind)    :: energy(nstate)  ! Internal energies
   NAMELIST/QNumbr/ mvib, jthrot, lorb, megaz, nuj, mtype ,energy
END MODULE Qnumber_Module
