MODULE Ngrid_Module
!  parameters for q1d to use in its search for 1D energy levels.

!  ndepth is the number of grid points with energy less than zero.
!  ngrid is the number of energy grid points with energy>zero.
!  mgrid is the maximum size of ngrid.  It sets dimensions.
!  scale1 determines the spacing of the grid for bound states.
!  scale2 determines the spacing of the grid for unbound states.      
   USE Numeric_Kinds_Module
   IMPLICIT NONE
   SAVE
   INTEGER(KIND=IW_Kind) ntheta !  ntheta   Number of theta angles.
   INTEGER(KIND=IW_Kind) nchi   !  nchi     Number of chi angles.
END MODULE Ngrid_Module
