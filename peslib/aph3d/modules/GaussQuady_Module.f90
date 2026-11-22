MODULE GaussQuady_Module
   USE Narran_Module
   USE Numeric_Kinds_Module
   USE Gaussb_Module
   SAVE
   REAL(Kind=WP_Kind) alpha(narran)              ! Harmonic oscillator exponential factor.
   REAL(Kind=WP_Kind), ALLOCATABLE:: xpth(:,:)   ! Hermite quadrature points.
   REAL(Kind=WP_Kind), ALLOCATABLE:: wpth(:,:)   ! Hermite quadrature weights.
END MODULE GaussQuady_Module
