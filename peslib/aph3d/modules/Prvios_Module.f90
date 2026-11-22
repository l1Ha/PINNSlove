MODULE Prvios_Module
   USE Numeric_Kinds_Module
   USE DVR_Module
   IMPLICIT NONE
   SAVE
   REAL(Kind=WP_Kind) tchip(nchimax,nchimax)      ! Chi-DVR transformation matrices for previous rho value
   REAL(Kind=WP_Kind) tthetap(nthetmax,nthetmax)  ! Theta-DVR transformation matrices for previous rho value
END MODULE Prvios_Module
