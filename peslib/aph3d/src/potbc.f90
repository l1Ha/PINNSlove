FUNCTION potbc (svib)
USE Narran_Module
USE FileUnits_Module
USE Arrch_Module
USE das_Module
USE Oops_Module
USE Masses_Module
USE potzero_Module
USE Numbers_Module
USE Numeric_Kinds_Module
IMPLICIT NONE
REAL(Kind=WP_Kind) cgama, potbc, rtrans, rvib, sj, svib, theta, v
INTEGER elect
!---------------------------------------------------------------------
! calculates the b-c potential.
!---------------------------------------------------------------------

sj=svib

IF(ioops)THEN
   theta=svib/rhoj
   sj=rhoj*sin(theta)
ENDIF

rvib=dscale(karran)*sj
rtrans=1.0d+30
cgama=0.0d0
CALL pothyp (rtrans,rvib,cgama,v)

IF(ioops)THEN
   v=v+0.25d0/(rhoj**2*usys2)
   IF(theta<0.0d0.or.theta>halfpi)THEN
      WRITE(Out_Unit,*)' warning rhoj,theta=',rhoj,theta
      v=500.0d0
   ENDIF
ENDIF

potbc=v

RETURN
END
