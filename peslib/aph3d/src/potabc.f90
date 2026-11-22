SUBROUTINE potabc (strans, svib, cgama, v, vp, vpp)
USE FileUnits_MODULE
USE Narran_Module
USE parms_MODULE
USE Oops_Module
USE Masses_Module
USE GaussQuady_Module
USE Arrch_Module
USE das_MODULE
USE Numbers_Module
USE Numbers_Module
IMPLICIT NONE
INTEGER index
REAL(Kind=WP_Kind) se, sbig, strans, slittl, svib, rho, theta
REAL(Kind=WP_Kind) rvib, rtrans, cgama, v, potbc, vp, vpp
!---------------------------------------------------------------------
! calculates the a-bc potential.
!---------------------------------------------------------------------

se=re(karran)
sbig=strans
slittl=svib

IF(ioops)THEN
    rho=strans
    theta=svib/strans
    sbig=rho*cos(theta)
    slittl=rho*sin(theta)
ENDIF

index=1
rvib=dscale(karran)*slittl
rtrans=sbig/dscale(karran)
CALL pothyp (rtrans,rvib,cgama,v)
v=v-potbc(svib)

IF(ioops)THEN
    v=v+quarter/(rho**2*usys2)
    IF(theta<zero.or.theta>halfpi)THEN
    WRITE(Out_Unit,*)' warning rho,theta=',rho,theta
    v=zero
    ENDIF
ENDIF

RETURN
!----------------***end-potabc***---------------------------------------
END
