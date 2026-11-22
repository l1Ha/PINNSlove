SUBROUTINE DVR_kinetic(ntheta,t_theta,tht_val,usys2)
!----------------------------------------------------------------------
! USE the DVR method to construct the theta kninetic energy
!---------------------------------------------------------------------
USE Numeric_Kinds_Module
USE FileUnits_Module
USE DVR2_MODULE

INTEGER :: ntheta, ltheta, i, j
REAL(Kind=WP_Kind) t_theta(ntheta,ntheta),tht_val(ntheta),usys2
LOGICAL debug
debug=.false.
ltheta=0
ntheta=ntheta-1
!
!     set up dvr transformation for theta
!
IF(Ntheta>NthetMax.or.Nchi>NchiMax)THEN
   WRITE(*,*)"Ntheta or Nchi exceeds maximum allowed values"
   WRITE(*,*)"Ntheta=",Ntheta,"  NthetMax=",NthetMax
   WRITE(*,*)"Nchi=",Nchi,"  NchiMax=",NchiMax
   STOP "Kinetic"
ENDIF
CALL tdv2fb1(ttheta, pttheta, whtheta, ntheta, ltheta, nthetmax, debug)
!
!     note that number of basis functions is ntheta-ltheta+1 for theta
!     Therefore, we readjust the value of
!     ntheta before we start the REAL calculation.
!
!
! DO this once and once only
!
ntheta=ntheta-ltheta+1
!
!
!     angmom determines the dvr j^2 matrix for theta and
!     puts it in antothet
!
CALL anmomthe(ttheta, pttheta, ntheta, nthetmax, ltheta, antothet, thetaj, debug)
!
!
antothet=antothet/usys2

!------------------
!get the theta values.
!------------------
DO i=1,ntheta
   tht_val(i)=acos(pttheta(i))/2.0d0
   DO j=1,ntheta
      t_theta(i,j)=antothet(i,j)
   ENDDO
ENDDO

!----------------
!check the sym
!----------------
DO i=1,ntheta
   DO j=1,i-1
      IF(ABS(t_theta(i,j)-t_theta(j,i))>1.0d-10)STOP 'wrong'
   ENDDO
ENDDO

RETURN
END
