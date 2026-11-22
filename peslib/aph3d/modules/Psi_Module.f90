MODULE Psi_Module
   USE Int1d_Module
   USE NState_Module
   IMPLICIT NONE
   SAVE
   REAL(Kind=WP_Kind) rho_basis                  !
   REAL(Kind=WP_Kind) eigen_old(nstate)          !
   REAL(Kind=WP_Kind) eigen_new(nstate)          !
   REAL(Kind=WP_Kind) psi_old(0:nltheta,nstate)  !
   REAL(Kind=WP_Kind) psi_new(0:nltheta,nstate)  !
   REAL(Kind=WP_Kind) Deau_new(Nstate)
   REAL(Kind=WP_Kind) D2eau_new(Nstate)

END MODULE Psi_Module
