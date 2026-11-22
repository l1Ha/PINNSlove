SUBROUTINE PotEff (rho, theta, xchi, potef)
USE Numeric_Kinds_Module
USE Parms_Module
USE FileUnits_Module
USE Masses_Module
USE TotJ_Module
USE MassFactor2_Module
USE Numbers_Module
!            P U R P O S E   O F    S U B R O U T I N E
!
! This routine calculates the effective potential (Hartree atomic
!  units).
!
!            I N P U T    A R G U M E N T S
!
!   rho   hyperspherical radius.              (0 < rho   <infinity)
!   theta hyperspherical polar angle theta.   (0 < theta < pi/2)
!   xchi  hyperspherical azimuthal angle chi. (0 < chi   < 2pi)
!
!            O U T P U T    A R G U M E N T S
!
!   potef  effective potential in Hartree atomic units
!          times 2*mu*rho*rho.

IMPLICIT NONE
INTEGER megas, ipfac, jjp
REAL(Kind=WP_Kind) pot, potef, rho,stheta,stsq, theta, xchi, tmurs
REAL(Kind=WP_Kind) ap, bp, cp, sin
EXTERNAL potaph
REAL(Kind=WP_Kind), PARAMETER:: eps=1.d-30

!-----------------------------------------------------------------------
! This CALL returns the interaction potential in Hartree atomic units.
!-----------------------------------------------------------------------
CALL potaph (rho, theta, xchi, pot)
!
stheta=sin(theta)
stsq=stheta*stheta
!-----------------------------------------------------------------------
! Evaluate the cs centrifugal potential.  For nonzero jtot and mega,
! this centrifugal potential is infinite at theta=pi/2 and theta=0.0.
! Hence, set the centrifugal potential equal to 1/eps there.
!  ap, bp, and cp are tmurs times the a, b, and c of the theory paper.
!-----------------------------------------------------------------------
ipfac=(-1)**(jtot+parity)
tmurs=usys2*rho*rho
jjp=jtot*(jtot+1)
megas=mega*mega
ap=one+stheta
bp=two*stsq
IF(bp<eps)bp=eps
cp=one-stheta
IF(cp<eps)cp=eps
ap=two/ap
bp=two/bp
cp=two/cp
potef=fifteen/four+half*(ap+bp)*(jjp-megas)+cp*megas
!----------------------------------------------------------------------
!  the following tests the adequacy of the cs basis for mega=1 and
!  cstest=.true.
!  this adds -(-1)**(j+p)*j(j+1)*(a-b)/4 to the usual cs centrifugal
!  potential.  The terms are collected to make the b singularity vanish
!  for jtot=mega=1 and parity=0 (even).
!-----------------------------------------------------------------------
IF(mega==1.AND.cstest)THEN
   potef=fifteen/four+cp+quarter*ap*((2-ipfac)*jjp-2)
   potef=potef+quarter*bp*((2+ipfac)*jjp-2)
ENDIF
!-----------------------------------------------------------------------
! The effective potential is the sum of the interaction potential and
! the centrifugal potential.
!-----------------------------------------------------------------------
potef=potef+tmurs*pot
!-----------------------------------------------------------------------
! The variable potef now contains 2*mu*rho**2*(the effective potential),
! all in Hartree atomic units.
!-----------------------------------------------------------------------
RETURN
ENDSUBROUTINE PotEff
