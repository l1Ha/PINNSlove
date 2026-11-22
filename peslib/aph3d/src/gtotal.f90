SUBROUTINE gtotal(theta, phi, gtot, gtheta, gphi)
USE FileUnits_Module
USE Narran_Module
USE MassFactor2_Module
USE opts_Module
!USE weight_Module
USE STST_Module
USE rhosur_Module

!----------------------------------------------------------------------
!  This routine calculates the scaling function.
!----------------------------------------------------------------------

!  I N P U T    A R G U M E N T S
!  theta
!  phi
!  gtot
!  gtheta
!  gphi


IMPLICIT NONE
LOGICAL full
INTEGER karran, k
REAL(Kind=WP_Kind) theta, phi, g, gtot, gth, gtheta, gph, gphi, stheta, sphi, denom
REAL(Kind=WP_Kind)   zero, one, two, half, sin, cos
PARAMETER (zero=0.d0, one=1.d0, two=2.d0, half=0.5d0)
DIMENSION stheta(narran), sphi(narran), g(narran), gth(narran), gph(narran)
EXTERNAL poteff

DATA full/.false./

IF(full)THEN
   WRITE(Out_Unit,*)'betaw=  ', betaw
   WRITE(Out_Unit,*)'seq=   ', seq
   WRITE(Out_Unit,*)'chncof=', chncof
   WRITE(Out_Unit,*)'chij',chij
   WRITE(Out_Unit,*)'rhosurf',rhosurf
ENDIF
!----------------------------------------------------------------------
!  Determine the interatomic distances s(karran)
!  which are stored in common save.
!----------------------------------------------------------------------
chi3(1)=phi
chi3(2)=chi3(1)-chij(2,1)
chi3(3)=chi3(1)-chij(3,1)
DO k=1,narran
   s(k)=rhosurf*sqrt(half-half*sin(theta)*cos(two*chi3(k)))
ENDDO
!----------------------------------------------------------------------
!  Calculating the weighting function.
!----------------------------------------------------------------------
gtot   = 0
gtheta = 0
gphi   = 0
DO karran = 1, narran
   denom          = (one-sin(theta)*cos(two*(chi3(karran))))**half
   IF(denom == zero) denom = 1.d-70
   stheta(karran) = -half*cos(theta)*cos(two*(chi3(karran)))*rhosurf/sqrt(two)/denom
   sphi(karran)   = sin(theta)*sin(two*(chi3(karran)))*rhosurf/sqrt(two)/denom
   IF(scheme == 1)THEN
      g(karran)   = exp(-betaw(karran)*(s(karran)-seq(karran))**2)*chncof(karran)
      gth(karran) = -two*stheta(karran)*g(karran)*betaw(karran)*(s(karran)-seq(karran))
      gph(karran) = -two*sphi(karran)*g(karran)*betaw(karran)*(s(karran)-seq(karran))
   ELSE
      g(karran)   = exp(-betaw(karran)*s(karran))*chncof(karran)
      gth(karran) = -one*stheta(karran)*g(karran)*betaw(karran)
      gph(karran) = -one*sphi(karran)*g(karran)*betaw(karran)
   ENDIF
   gtot        = gtot + g(karran)
   gtheta      = gtheta + gth(karran)
   gphi        = gphi + gph(karran)
ENDDO
IF(gtot == zero)gtot=1.d-90

IF(full)THEN
   WRITE(Out_Unit,*)' rhosurf=', rhosurf, ' theta=', theta, ' phi=', phi
   WRITE(Out_Unit,*)'s=', s
   WRITE(Out_Unit,*)'stheta=', stheta
   WRITE(Out_Unit,*)'sphi=', sphi
   WRITE(Out_Unit,*)'seq=', seq
   WRITE(Out_Unit,*)'g=', g
   WRITE(Out_Unit,*)'gth=', gth
   WRITE(Out_Unit,*)'gph=', gph
   WRITE(Out_Unit,*)'gtot=', gtot, ' gtheta=',  gtheta, ' gphi=', gphi
   STOP 'gtotal'
ENDIF
 7    FORMAT(1x, 5e14.7)
RETURN
END
