SUBROUTINE Main1D(rho,potmax,nchi,ltheta,nch,theta, peigv,itheta,debug,h3sys,ii)
USE Numbers_Module
USE FileUnits_Module
USE Convrsns_Module
USE DVR_Module
USE DVRMasses_Module
USE OneD_Module
USE TotJ_Module
USE DVR2_Module
USE CSBasis_Module
USE Parms_Module
!
!     calculate 1-d eigenvalues and eigenvectors for constant theta
!     and rho.
!
IMPLICIT NONE
INTEGER nvmax, itheta, nctem,ichi, nchi, ii,i, j, jchi, ipfac, jjp
INTEGER k, ltheta, nch, megas
LOGICAL peigv, debug, h3sys
REAL(Kind=WP_Kind) dist(3)
REAL(Kind=WP_Kind) sinth, theta, rho, r, d, chi1, Potntl_Min
REAL(Kind=WP_Kind) potntl, pottem, potmax, pottemp, chii, CHI1J, temp1
REAL(Kind=WP_Kind) rtwomu, chi, eps, stheta, ap, bp, cp, poteff, sd, subd, en
!
!     a few temporary matrices. subd and sd are used to diagonalize
!     the hamiltonian. pottemp keeps track of the values of the
!     potential energy at the dvr points
!
COMMON /waste/subd(nmax),  sd(nmax*(nmax+1)/2),  en(nmax),pottemp(nchimax)
!
nvmax=nmax*(nmax+1)/2
nev(itheta)=0
sinth=1.d0/sin(theta/2)**2
!
!     nctem will be the number of points kept after truncation
!
nctem=1
rtwomu=1.d0/(twomu*rho*rho)
!
!     loop over all the chi dvr points and calculate the potential
!     at each point. IF the energy is less than potmax, include
!     the point in the basis
!
IF(Nchi>NchiMax)THEN
   WRITE(*,*)"Ntheta or Nchi exceeds maximum allowed values"
   WRITE(*,*)"Nchi=",Nchi,"  NchiMax=",NchiMax
   STOP "Main1D"
ENDIF
DO 1 ichi=1,nchi
   IF(nctem>nmax)THEN
      WRITE(Out_Unit,'(1x," nctot = ",i4)')nctem
      STOP
   ENDIF
   chi=acos(ptchi(ichi))
   CALL pptorr(rho,theta/Two,chi,r,d,chi1)
   CALL rds2r3(r,d,chi1,dist,nch)
!
   IF(ii==1)THEN
     CALL surface(potntl,dist)
   ENDIF
!
   pottem=potntl*autoev
!
!     because pk2 surface has some negative values when atoms get very
!     CLOSE together, we have extra condition on pottem
!
   IF( (h3sys) .AND.( (pottem>potmax).or.(pottem<=0.0d0)) ) goto 1
!
!     the above condition applicable only to h3 pk2 surface which
!     has some negative eigenvalues when atoms are very CLOSE. in other
!     cases (surfaces) the condition below applies
!
  IF( (.NOT.h3sys) .AND. (pottem>potmax) ) goto 1
   pottemp(nctem)=pottem
!
!     nptchi keeps track tells us that basis point nctem corresponds
!     to dvr point ichi at this value of theta
!
   nptchi(nctem,itheta)=ichi
   nctem=nctem+1
 1    CONTINUE
!
!     the following line is needed to get the number of angles right
!
nctem=nctem-1
nc2tot(itheta)=nctem
!
!     now we calculate the hamiltonian by the usual methods
!
Potntl_Min=HUGE(One)
DO i=1,nctem
   ichi=nptchi(i,itheta)
   chii=acos(ptchi(ichi))
   DO j=i,nctem
      jchi=nptchi(j,itheta)
      chi1j=acos(ptchi(jchi))
      temp1=0.0d0
      !
      !     calculate potential only for the diagonal matrix elements
      !
      IF(i==j)THEN
      !
      !     calculate the interaction potential
      !
         CALL pptorr(rho,theta/Two,chii,r,d,chi1)
         CALL rds2r3(r,d,chi1,dist,nch)
      !
         IF(ii==1)THEN
           CALL surface(potntl,dist)
         ENDIF
      !
      !
      !     calculate the cs centrifugal potential
      !
         IF(jtot/=0)THEN
           ipfac=(-1)**(jtot+parity)
           jjp=jtot*(jtot+1)
           megas=mega*mega
           eps=1.0d-30
           stheta=sin(0.5d0*theta)
           ap=1.0d0+stheta
           bp=2.0d0*stheta*stheta
           IF(bp<eps)bp=eps
           cp=1.0d0-stheta
           IF(cp<eps)cp=eps
           ap=2.0d0/ap
           bp=2.0d0/bp
           cp=2.0d0/cp
           IF(.NOT.cstest)THEN
             poteff=0.5d0*(ap+bp)*(jjp-megas) + cp*megas
           ELSE
             poteff=cp + 0.25d0*ap*((2-ipfac)*jjp-2)+ 0.25d0*bp*((2+ipfac)*jjp-2)
           ENDIF
           temp1=potntl*autoev + poteff*rtwomu
         ELSE
           temp1=potntl*autoev
         ENDIF
      ENDIF
      !
      !     now add the dvr kinetic energy element multiplied by
      !     1./sin(theta/2)**2 to temp1 to get the matrix element
      !
      h1d(i,j,itheta)=rtwomu*((antotchi(ichi,jchi))*sinth)+temp1
      h1d(j,i,itheta)=h1d(i,j,itheta)
      IF(Potntl_Min>Temp1)Potntl_Min=Temp1
   ENDDO
ENDDO
!WRITE(*,*)"Potntl_Min=",Potntl_Min  !TMPMODGregParker needs to be passed to eliminate some eigenvectors
!
!     stuff the lower triangle of h1d into sd and CALL tdiag to give
!     eigenvalues and eigenvectors
!
k=0
DO i=1,nctem
   DO j=1,i
      k=k+1
      sd(k)=h1d(i,j,itheta)
   ENDDO
ENDDO
WRITE(Out_Unit,*)'nctem,nmax,nvmax=',nctem,nmax,nvmax
IF(nctem>=1.)THEN
   CALL tdiagrw(sd,e1d(1,itheta),cc1dtot(1,1,itheta),subd,nctem,nmax,nvmax)
ENDIF
!
RETURN
ENDSUBROUTINE Main1D
