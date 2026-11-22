MODULE Contour_Module
   USE Numeric_Kinds_Module
   IMPLICIT NONE
   SAVE
   INTEGER(KIND=IW_Kind) :: nx  =301    ! Number of x grid points
   INTEGER(KIND=IW_Kind) :: ny  =301    ! number of y grid points
   INTEGER(KIND=IW_Kind) :: iaph=1      ! Surface function index
   REAL(Kind=WP_Kind)    :: rho =1.d0   ! hyperradius
   REAL(Kind=WP_Kind)    :: VMax=4.0D0  ! Maximum Potential used in Plotting
   NAMELIST/contour/ nx, ny, iaph, rho, VMax
END MODULE Contour_Module
