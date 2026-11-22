SUBROUTINE grms(theta, phi, gtot, gtheta, gphi)
USE FileUnits_Module
USE MassFactor2_Module
USE opts_Module
USE STST_Module
USE rhosur_Module

!  P U R P O S E   O F    S U B R O U T I N E

!----------------------------------------------------------------------
!  This routine calculates the scaling function.
!----------------------------------------------------------------------

!  I N P U T    A R G U M E N T S
!  theta
!  phi
!  gtot
!  gtheta
!  gphi



!  L O G I C A L S

LOGICAL full

!  I N T E G E R S

INTEGER  k

!  R E A L S
REAL(Kind=WP_Kind) theta, phi, g, gtot, gth, gtheta, gph, gphi, stheta, sphi, denom
REAL(Kind=WP_Kind)   zero, one, two, half
PARAMETER (zero=0.d0, one=1.d0, two=2.d0, half=0.5d0)

!  D I M E N S I O N S

DIMENSION stheta(narran), sphi(narran), g(narran), gth(narran), gph(narran)

!  I N T R I N S I C    F U N C T I O N S
INTRINSIC sin, cos

!  E X T E R N A L S
EXTERNAL poteff

DATA full/.false./

IF(full)THEN
   WRITE(Out_Unit,*)'betaw=  ', betaw
   WRITE(Out_Unit,*)'seq=   ', seq
   WRITE(Out_Unit,*)'chncof=', chncof
ENDIF
!----------------------------------------------------------------------
!  Determine the interatomic distances s(karran)
!  which are stored in common save.
!----------------------------------------------------------------------
chi3(1)=phi
chi3(2)=chi3(1)-chij(2,1)
chi3(3)=chi3(1)-chij(3,1)
!----------------------------------------------------------------------
!  Calculating the weighting function.
!----------------------------------------------------------------------
gtot   = 0
gtheta = 0
gphi   = 0
DO 1 k = 1, narran
   s(k)=rhosurf*sqrt(half-half*sin(theta)*cos(two*chi3(k)))
   IF(scheme == 1)THEN
      g(k)   = exp(-betaw(k)*(s(k)-seq(k))**2)*chncof(k)
   ELSE
      g(k)   = exp(-betaw(k)*s(k))*chncof(k)
   ENDIF
   gtot        = gtot + g(k)
 1    CONTINUE
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
   WRITE(Out_Unit,*)'gtot=', gtot, ' gtheta=', gtheta, ' gphi=', gphi
ENDIF
 7    FORMAT(1x, 5e14.7)
RETURN
END
