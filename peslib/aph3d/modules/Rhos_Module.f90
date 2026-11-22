MODULE Rhos_Module
   USE Numeric_Kinds_Module
   IMPLICIT NONE
   SAVE
   REAL(Kind=WP_Kind) rho       ! current rho value.
   REAL(Kind=WP_Kind) rhomin    ! Minimum rho value.
   REAL(Kind=WP_Kind) rhomax    ! Maximum rho value.
   REAL(Kind=WP_Kind) deltarho  ! Rho spacing.
END MODULE Rhos_Module
