MODULE DVR2_Module
   USE Numeric_Kinds_Module
   USE DVR_Module
   IMPLICIT NONE
   SAVE
   REAL(Kind=WP_Kind) tchi(nchimax,nchimax)       ! chi-dvr transformation matrices
   REAL(Kind=WP_Kind) ttheta(nthetmax,nthetmax)   ! theta-dvr transformation matrices
   REAL(Kind=WP_Kind) ptchi(nchimax)              ! chi-dvr points
   REAL(Kind=WP_Kind) pttheta(nthetmax)           ! theta--dvr points
   REAL(Kind=WP_Kind) whchi(nchimax)              ! whchi-here for historical reasons
   REAL(Kind=WP_Kind) whtheta(nthetmax)           ! whtheta--here for historical reasons
   REAL(Kind=WP_Kind) chij(nchimax)               ! temporary storage for calculation of dvr j^2 matrices
   REAL(Kind=WP_Kind) chij1(nchimax,nchimax)      ! temporary storage for calculation of dvr j^2 matrices
   REAL(Kind=WP_Kind) thetaj(nthetmax)            ! temporary storage for calculation of dvr j^2 matrices
   REAL(Kind=WP_Kind) antotchi(nchimax,nchimax)   ! antotchi-dvr j^2 matrices
   REAL(Kind=WP_Kind) antothet(nthetmax,nthetmax) ! antothet--dvr j^2 matrices
END MODULE DVR2_Module
