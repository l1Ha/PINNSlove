SUBROUTINE rrtopp(r,d,chi,rho,theta,phi)
USE Numeric_Kinds_Module
!
! transforms Jacobi coordinates r,d,chi to Pack's and Parker's hyperspherical
! coordinates rho,theta,phi. reference:Pack and Parker reprint
!
REAL(Kind=WP_Kind) rho, r, d, temp, chi, theta, phi, rs2, st, cp
rho=sqrt(r*r+d*d)
temp=(r*r-d*d)**2+(2.d0*r*d*cos(chi))**2
theta=atan2(sqrt(temp),(2.d0*r*d*sin(chi)))
phi=0.5d0*atan2(2.d0*r*d*cos(chi),r*r-d*d)
RETURN
!
entry pptorr(rho,theta,phi,r,d,chi)
!
!     inverse transformation between Pack's and Parker's coordinates and Jacobi
!     coordinates
!
rs2=1.d0/sqrt(2.d0)
st=sin(theta)
cp=cos(2*phi)
r=rho*rs2*sqrt(1.d0+st*cp)
d=rho*rs2*sqrt(1.d0-st*cp)
chi=acos(st*sin(2.d0*phi)/sqrt(1.d0-st*st*cp*cp))
RETURN
!
END
