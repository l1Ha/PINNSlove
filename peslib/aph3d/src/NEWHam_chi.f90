SUBROUTINE NewHam_chi(t_chi,nchir,itheta,B_theta,rhom2,potmax,v_pot,ntheta,nrho,irho)
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
USE TwoD_Module,h2d=>csurf1

!
!     calculate 1-d eigenvalues and eigenvectors for constant theta
!     and rho.
!
IMPLICIT NONE
INTEGER nrho, ntheta, nchir, nvmax, itheta, nctem, ichi, irho, i, j, jchi, k
LOGICAL h3sys
REAL(Kind=WP_Kind) v_pot(nrho,ntheta,nchir), t_chi(nchir,nchir)
REAL(Kind=WP_Kind) pottem, potmax, temp1, b_theta, rhom2
REAL(Kind=WP_Kind) sd, subd, en, pottemp
!
!     h2d--2-d hamiltonian
!
!COMMON/twod/ h2d(nt2max,nt2max),       csurf2(ndvptmax,nsfnmax), ovlp(nsfnmax,nsfnmax)
!
!
!     a few temporary matrices. subd and sd are used to diagonalize
!     the hamiltonian. pottemp keeps track of the values of the
!     potential energy at the dvr points
!
COMMON/waste/subd(nmax),  sd(nmax*(nmax+1)/2),  en(nmax), pottemp(nchimax)




!
h3sys=.false.
nvmax=nmax*(nmax+1)/2
nev(itheta)=0

!
!     nctem will be the number of points kept after truncation
!
nctem=1
!
!     loop over all the chi dvr points and calculate the potential
!     at each point. IF the energy is less than potmax, include
!     the point in the basis
!
DO 1 ichi=1,nchir
   IF(nctem>nmax)THEN
      WRITE(Out_Unit,'(1x," nctot = ",i4)')nctem
      STOP "NewHam_Chi"
   ENDIF
   pottem=V_pot(irho,itheta,ichi)*autoev
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
DO i=1,nctem
   ichi=nptchi(i,itheta)
   DO j=i,nctem
      jchi=nptchi(j,itheta)
      temp1=0.0d0
      !
      !     calculate potential only for the diagonal matrix elements
      !
      IF(i==j)THEN
            temp1=v_pot(irho,itheta,ichi)
      ENDIF
      !
      !     now add the dvr kinetic energy element multiplied by
      !     1./sin(theta/2)**2 to temp1 to get the matrix element
      !
      h1d(i,j,itheta)=t_chi(ichi,jchi)*2.d0 * b_theta*rhom2+temp1
      h1d(j,i,itheta)=h1d(i,j,itheta)
   ENDDO
ENDDO
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

RETURN
END
