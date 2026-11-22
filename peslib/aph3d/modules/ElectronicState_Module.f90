MODULE ElectronicState_Module
   USE Numeric_Kinds_Module
   IMPLICIT NONE
   SAVE
   INTEGER :: IthElec, KthElec
   REAL(KIND=WP_Kind) :: hstep
   REAL(KIND=WP_Kind), ALLOCATABLE:: DiatomicPotential(:,:), Dipole(:,:,:), rvals(:)
END MODULE ElectronicState_Module