MODULE Ngrid4_Module
!  parameters for q1d to use in its search for 1D energy levels.
!  ndepth is the number of grid points with energy less than zero.
!  ngrid is the number of energy grid points with energy>zero.
!  mgrid is the maximum size of ngrid.  It sets dimensions.
!  scale1 determines the spacing of the grid for bound states.
!  scale2 determines the spacing of the grid for unbound states.      
   USE Numeric_Kinds_Module
   IMPLICIT NONE
   SAVE
   INTEGER(KIND=IW_Kind), PARAMETER:: NGrid=600
   INTEGER(KIND=IW_Kind), PARAMETER:: MGrid=4000
   INTEGER(KIND=IW_Kind), PARAMETER:: NDepth=-120
   REAL(Kind=WP_Kind), PARAMETER:: Scale1=0.01414141414141414d0
   REAL(Kind=WP_Kind), PARAMETER:: Scale2=0.05d0
END MODULE Ngrid4_Module
