MODULE Int1D_Module
   USE Numeric_Kinds_Module
   IMPLICIT NONE
   SAVE
   INTEGER(KIND=IW_Kind), PARAMETER:: nbtheta=100    !
   INTEGER(KIND=IW_Kind), PARAMETER:: mxjrot=16      !
   INTEGER(KIND=IW_Kind), PARAMETER:: n1=700         !
   INTEGER(KIND=IW_Kind), PARAMETER:: n2=1000        !
   INTEGER(KIND=IW_Kind), PARAMETER:: nltheta=n1+n2   !
   INTEGER(KIND=IW_Kind), PARAMETER:: itmax=100      !
   REAL(Kind=WP_Kind), PARAMETER:: r0=3.1d0       !
   REAL(Kind=WP_Kind), PARAMETER:: rswitch=17.0d0 !
   REAL(Kind=WP_Kind), PARAMETER:: tol=1.0d-13    !
   REAL(Kind=WP_Kind) switch                      !
   REAL(Kind=WP_Kind) lthetzro                    !
   REAL(Kind=WP_Kind) sw_old                      !
   REAL(Kind=WP_Kind) ltz_old                     !
END MODULE Int1D_Module
