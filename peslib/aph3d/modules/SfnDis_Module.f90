MODULE SfnDis_Module
   USE Numeric_Kinds_Module
   IMPLICIT NONE
   SAVE
   INTEGER(KIND=IW_Kind) :: ithrho     !
   INTEGER(KIND=IW_Kind) :: ithmax     !
   INTEGER(KIND=IW_Kind) :: NRho      = 301       !
   REAL(Kind=WP_Kind)    :: rhomin    = 2.1d0     ! First rho value for calculating surface functions.
   REAL(Kind=WP_Kind)    :: rhomax    = 8.10D0    ! Last rho value for calculating surface functions.
   REAL(Kind=WP_Kind)    :: deltarho1 = 0.05d0  !
   REAL(Kind=WP_Kind)    :: deltarho2 = 0.d0  !
   REAL(Kind=WP_Kind)    :: rho2nd     !
   REAL(Kind=WP_Kind)    :: rhostop    !
   NAMELIST/sfundist/ deltarho1, deltarho2, rhomin, rhomax, nrho
END MODULE SfnDis_Module
