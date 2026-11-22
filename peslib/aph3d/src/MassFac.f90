SUBROUTINE MassFac(little)
!-----------------------------------------------------------------------
!  Routine to calculate all the factors that are only dependent
!  on the masses.
!-----------------------------------------------------------------------
USE Numeric_Kinds_Module
USE Parms_Module
USE FileUnits_Module
USE Masses_Module
USE MassFactor2_Module
USE Convrsns_Module
USE ChiAng_Module
USE Das_Module
USE Numbers_Module
IMPLICIT NONE
REAL(Kind=WP_Kind) totmas, cchi, schi, chiab, chibc, chica, angp
REAL(Kind=WP_Kind) angm
INTEGER i, ip, im
LOGICAL little
!-----------------------------------------------------------------------
mass(1)=mass(1)*amutoau
mass(2)=mass(2)*amutoau
mass(3)=mass(3)*amutoau
totmas=mass(1)+mass(2)+mass(3)
UBC=mass(2)*mass(3)/(mass(2)+mass(3))
UAB=mass(1)*mass(2)/(mass(1)+mass(2))
UAC=mass(1)*mass(3)/(mass(1)+mass(3))
U_Diatom(1)=UBC
U_Diatom(2)=UAC
U_Diatom(3)=UAB
usys=SQRT(mass(1)*mass(2)*mass(3)/totmas)
usys2=two*usys
UABC=usys
UABC2=usys2
mtot=totmas
dscale(1)=SQRT(mass(1)/usys*(one-(mass(1)/totmas)))
dscale(2)=SQRT(mass(2)/usys*(one-(mass(2)/totmas)))
dscale(3)=SQRT(mass(3)/usys*(one-(mass(3)/totmas)))
cchi=-usys/(dscale(1)*dscale(2)*mass(3))
schi=one/(dscale(1)*dscale(2))
chiab=ATAN2(schi,cchi)
cchi=-usys/(dscale(2)*dscale(3)*mass(1))
schi=one/(dscale(2)*dscale(3))
chibc=ATAN2(schi,cchi)
cchi=-usys/(dscale(3)*dscale(1)*mass(2))
schi=one/(dscale(3)*dscale(1))
chica=ATAN2(schi,cchi)
cchi=chiab+chibc+chica
chif1(1)=Zero
chif1(2)=-chiab
chif1(3)=chica
IF(little)THEN
   WRITE(Out_Unit,*)
   WRITE(Out_Unit,*)'MassFac: Start'
   WRITE(Out_Unit,310) mass(1), mass(2), mass(3)
   WRITE(Out_Unit,311) totmas, usys
   WRITE(Out_Unit,312) uabc, uabc2, ubc
   WRITE(Out_Unit,320) chiab, chibc, chica, cchi
   WRITE(Out_Unit,321) chif1
   WRITE(Out_Unit,322) dscale
ENDIF
312 FORMAT(5x,'uabc =',ES14.7,' uabc2=',ES14.7, '   ubc=',ES14.7)
310 FORMAT(5x,'mass(1)=',ES14.7,' mass(2)=', ES14.7,' mass(3)=',ES14.7)
311 FORMAT(5x,'Total mass=', ES14.7,' Hyperspherical reduced mass=', ES14.7)
320 FORMAT(5x,'chi-ab=', ES14.7, ' chi-bc=', ES14.7,' chi-ca=', ES14.7, ' chi-sum=', ES14.7)
321 FORMAT(5x,'chif1 =',ES14.7)
322 FORMAT(5x,'dscale=',ES14.7)
!-----------------------------------------------------------------------
!  More mass factors.
!  chij(i+1,i) is negative obtuse, kinematic angle between channels.
!-----------------------------------------------------------------------
xmass(1)=mass(1)
xmass(2)=mass(2)
xmass(3)=mass(3)
DO i=1, narran
   chij(i,i)=Zero
   d(i)=dscale(i)
ENDDO
DO i=1, narran
   ip=i+1
   im=i-1
   IF(i == narran)ip=1
   IF(i == 1)im=narran
   chij(i,ip)=ACOS(-usys/(dscale(i)*dscale(ip)*xmass(im)))
   chij(ip,i)=-chij(i,ip)
ENDDO
!  ------------------------------------------------------------
!  acute(i) is the smallest of the acute angles between channel i
!  and its neighbors.  jacute(i) is the index of that nearest neighbor.
!  --------------------------------------------------------------------
DO i=1, narran
   ip=i+1
   im=i-1
   IF(i == narran)ip=1
   IF(i == 1)im=narran
   angp=chij(ip,i)+pi
   angm=chij(i,im)+pi
   IF(angp <= angm)THEN
      acute(i)=angp
      jacute(i)=ip
   ELSE
      acute(i)=angm
      jacute(i)=im
   ENDIF
ENDDO

IF(Little)THEN
   WRITE(Out_Unit,330) (jacute(i), acute(i), i=1, narran)
   WRITE(Out_Unit,*)'MassFac: End'
   WRITE(Out_Unit,*)
ENDIF
!-----------***END-massfac***-------------------------------------------
RETURN
330 FORMAT(1x, 'jacute, acute=', 3(i5, 1x, es14.7))
ENDSUBROUTINE MassFac
