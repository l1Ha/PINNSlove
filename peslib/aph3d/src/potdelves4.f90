REAL(Kind=WP_Kind) FUNCTION potdelves4(rho, ltheta, btheta)
!-------------------------------------------------------------------
! funcrtion : potdelves4
!
! package   : CID
!
! Language  : Fortran 90
!
! author    :  R. T Pack
!           F. Colavecchia (flavioc@lanl.gov)
!
! date      :              version:
! revision  :              version:
!
! purpose   :   This routine generates the perturbing potential needed in
!               the Delves CC equations.  See RTP notes of 98/4/20-3 Eq. (20).
!               By specifying various functional forms and doing the integral
!               analytically one can compare with the calculated numerical
!               results for accuracy and rates of convergence.
!
!
! input     :   iarr -> arrangement number
!               rho  -> \
!              ltheta-> - Hyperspherical coordinates.
!              btheta-> /
! output    :      ->
!                  ->
!
! modules   :
!
!
! common    :
!
!
! notes     :
!
!-------------------------------------------------------------------
USE Masses_Module
USE Das_Module
USE psi_Module
USE Narran_Module
USE int1d_Module
IMPLICIT NONE
CHARACTER(LEN=21), PARAMETER:: ProcName='potdelves4'
REAL(Kind=WP_Kind) rho, ltheta, btheta
REAL(Kind=WP_Kind) ss, rrr, cosf, xt, xm(3), rij(3), s(3), vpot, v3pot

! This routine generates the perturbing potential needed in
! the Delves CC equations.  See RTP notes of 98/4/20-3 Eq. (20).
! By specifying various functional forms and doing the integral
! analytically one can compare with the calculated numerical
! results for accuracy and rates of convergence.

potdelves4=0.d0
IF(potdelves4/=-123.d0)THEN
   STOP "THIS ROUTINE IS NOT WOEKING"
ENDIF
xm(1)=mass(1)
xm(2)=mass(2)
xm(3)=mass(3)
!  get scaled and unscaled interparticle distances.
!  distances of channel 1
!  these go with the ACTUAL rho.
s(1) = rho*sin(ltheta)
ss = rho*cos(ltheta)
rij(1) = dscale(1)*s(1)
rrr = ss/dscale(1)
!  distances of other channels
cosf = cos(btheta)
!  temporary quantity
xt = rij(1)*xm(2)/(xm(2)+xm(3))
rij(2) = sqrt(rrr**2+xt**2-2.d0*rrr*xt*cosf)
xt = rij(1)*xm(3)/(xm(2)+xm(3))
rij(3) = sqrt(rrr**2+xt**2+2.d0*rrr*xt*cosf)
!CALL surfaceneh2(vpot, rij)
!  The above is the complete PES
!  now get the v_3 for SECTOR CENTER rho and
!  subtract it off.
s(1) = rho_basis*sin(ltheta)
rij(1) = dscale(1)*s(1)
!  surfaca only uses rij(1), so we can skip the following several
!  statements
!      ss = rho_basis*cos(ltheta)
!      rrr = ss/dscale(1)
!  distances of other channels
!      cosf = cos(btheta)
!  temporary quantity
!      xt = rij(1)*xm(2)/(xm(2)+xm(3))
!      rij(2) = sqrt(rrr**2+xt**2-2.d0*rrr*xt*cosf)
!      xt = rij(1)*xm(3)/(xm(2)+xm(3))
!      rij(3) = sqrt(rrr**2+xt**2+2.d0*rrr*xt*cosf)
!CALL surfaca(v3pot, rij)
!
potdelves4=vpot-(rho_basis/rho)**2*v3pot
! The above statement returns the integrand of RTP 98/4/30-3
! equation 20, not the interaction potential.

RETURN
END FUNCTION potdelves4
