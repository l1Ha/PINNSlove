SUBROUTINE NEWHam_2d(rhom2,ntheta,nchi,t_theta,irho,lam,jtot,a,b,c,parity,nbig2)
USE FileUnits_Module
USE DVR_Module
USE DVRMasses_Module
USE OneD_Module
USE DVR2_Module
USE Parms_Module
USE TwoD_Module,h2d=>csurf1
USE Masses_Module
!
!     this SUBROUTINE takes the 1-d hamiltonians from hmat1d and the
!     coupling matrices from angmom and turns them into a full 2-d
!     hamiltonian matrix.
!     the resulting matrix is stored in h2d which
!     is THEN diagonalized to give eigenvectors (in h2d) and
!     eigenvalues (in e2d).
!
IMPLICIT NONE
INTEGER ntheta, ntvmax, ngri, ltheta, nchi, itheta, i, j, ichi, ngrj, jchi, indexi
INTEGER ntempi, jtheta, ntempj, jtot, lam, i4, item4, i1, i2, k, in, item1, itemp
REAL(Kind=WP_Kind) t_theta(ntheta,ntheta),rhom2,a(ntheta),b(ntheta),c(ntheta)
REAL(Kind=WP_Kind),ALLOCATABLE ::temparray(:,:),temparray2(:,:)
REAL(Kind=WP_Kind) rtwomu, factor, cfactor, cfac, factor2, waa, htemp, rho, eff
REAL(Kind=WP_Kind) sd, subd
INTEGER irho, parity, nbig2
EXTERNAL mprint
PARAMETER (ntvmax= nt2max*(nt2max+1)/2)


!
!     h2d--2-d hamiltonian
!     later, when h2d is no longer needed, the final 2-d eigenvectors
!     are stored in h2d.
!
!COMMON/twod/ h2d(nt2max,nt2max),       csurf2(ndvptmax,nsfnmax), ovlp(nsfnmax,nsfnmax)
!
COMMON/waste/sd(ntvmax),    subd(nt2max)
rtwomu=1./(twomu*rho*rho)
ngri=0
ltheta=0

ALLOCATE(h2d(ntheta*nchi,ntheta*nchi))

IF(Ntheta>NthetMax.or.Nchi>NchiMax)THEN
   WRITE(*,*)"Ntheta or Nchi exceeds maximum allowed values"
   WRITE(*,*)"Ntheta=",Ntheta,"  NthetMax=",NthetMax
   WRITE(*,*)"Nchi=",Nchi,"  NchiMax=",NchiMax
   STOP "NewHam2D"
ENDIF

!
!     cc1dtot is in the truncated dvr basis. in order to make
!     the fast transformation work, we need to have a eigenvector
!     matrix in the full dvr basis. the following loop does that
!     putting the full eigenvector matrix in cc1d
!
DO itheta=1,ntheta
   DO i=1,nchimax
      DO j=1,nchimax
         cc1d(i,j,itheta)=0.0d0
      ENDDO
   ENDDO
ENDDO
!

DO itheta=1,ntheta
   DO i=1,nc2tot(itheta)
      ichi=nptchi(i,itheta)
      DO j=1,nc2tot(itheta)
         cc1d(ichi,j,itheta)=cc1dtot(i,j,itheta)
      ENDDO
   ENDDO
ENDDO
!
!     begin outer loop over theta
!

DO itheta=1,ntheta
   !
   !     ngrj and ngri will help us keep track of where we are in
   !     the full hamiltonian matrix. basically, they are partial
   !     sums of the number of 1-d eigenvectors
   !
   ngrj=ngri
   ntempi=nev(itheta)
   !
   !     begin inner loop over theta
   !
   DO jtheta=itheta,ntheta
      ntempj=nev(jtheta)
      eff=0.0d0
      !
      !     for diagonal elements (both theta points and chi points
      !     are the same), must add in effective potential
      !
      IF(itheta==jtheta)eff=3.75d0!-16d0*ltheta*ltheta*2.0d0*b(itheta)
      !-------------------------------------------------------------------
      !factor is the diagonal part when Jtot/=0
      !-------------------------------------------------------------------
      IF(itheta==jtheta)factor=jtot * (jtot + 1) * (a (itheta) + b (itheta) )  / usys2
      IF(itheta==jtheta)factor=factor+(2.d0 * c (itheta) - (a (itheta) + b (itheta) ) ) * lam**2 / usys2
      !-------------------------------------------------------------------
      ! construct Asymmetric Top coupling term for lambda=lambda'=1
      ! This is a diagonal contribution for only lambda=1.
      !-------------------------------------------------------------------
      cfactor = cfac (jtot, 1, 1, parity)
      IF(Jtot/=0..AND.itheta==jtheta.AND.lam==1)factor2=-(a (itheta) - b (itheta) ) / 2.d0 * cfactor/ usys2
      waa=(t_theta(itheta,jtheta)+eff/ usys2+factor+factor2)* rhom2
      !
      !     now transform the (diagonal) block of the coupling matrix
      !
      DO i4=1,ntempi
         item4=ngri+i4
         DO i1=1,ntempj
            item1=ngrj+i1
            htemp=0.0d0
            !
            !     because waa is diagonal in chi, we only need to DO one loop
            !     over chi. this is the heart of the fast transformation
            !
            DO i2=1,nchi
               htemp=htemp+cc1d(i2,i4,itheta)*cc1d(i2,i1,jtheta)*waa
            ENDDO
            !
            !     when done, place terms in the correct place in h2d
            !
            h2d(item4,item1)=htemp
            h2d(item1,item4)=htemp
         ENDDO
      ENDDO
      ngrj=ngrj+ntempj
   ENDDO
   !
   !     another time saving maneuver. don't add in full 1-d hamiltonian
   !     just add in 1-d eigenvalues along the diagonal
   !
   DO in=1,ntempi
      item1=ngri+in
      h2d(item1,item1)=h2d(item1,item1)+e1d(in,itheta)
   ENDDO
   ngri=ngri+ntempi
ENDDO
!
!     done constructing h2d
!
!
!     diagonalize h2d, the lower triangle of which is stored in the
!     waste array sd. after diagonalization, e2d will contain the
!     eigenvalues and the eigenvectors will be stored columnwise in
!     h2d.
!
nc3tot(irho)=ngri
ngrj=ngri
k=0
DO i=1,ngri
   DO j=1,i
      k=k+1
      sd(k)=h2d(i,j)
   ENDDO
ENDDO
h2d=0.0d0


ALLOCATE(temparray(ngri,ngri))

CALL tdiagrw(sd,e2d_new(1,irho),temparray,subd,ngri,ngri,ngri*(ngri+1)/2)

ALLOCATE(temparray2(ntheta*nchi,ntheta*nchi))

DO i=1,ngri
   DO j=1,ngri
      temparray2(i,j)=temparray(i,j)
   ENDDO
ENDDO

DEALLOCATE(temparray)

DO i=1,ngri
   DO j=1,ngri
      htemp=0.0d0
      DO k=1,ngri
      htemp=htemp+temparray2(k,i)*temparray2(k,j)
      ENDDO
      IF(i==j)THEN
           IF(ABS(htemp-1.0d0)>1.0d-6)THEN
           WRITE(Out_Unit, * )'wrong 2dbasis=',i,j,htemp
           STOP 'error 2d'
           ENDIF
      ELSE
           IF(ABS(htemp)>1.0d-6)THEN
           WRITE(Out_Unit, * )'wrong 2dbasis=',i,j,htemp
           STOP 'error 2d'
           ENDIF
      ENDIF
   ENDDO
ENDDO

DEALLOCATE(h2d)


!--------------------------------------
!NEW transformation matrix=d*C
!--------------------------------------
ALLOCATE(temparray(ntheta*nchi,ntheta*nchi))
DO k=1,ngri
   itemp=0
   DO itheta=1,ntheta
      DO ichi=1,nchi
         indexi=nchi*(itheta-1)+ichi
         DO jchi=1,nev(itheta)
            temparray(indexi,k)=temparray(indexi,k)+cc1d(ichi,jchi,itheta)*temparray2(itemp+jchi,k)
         ENDDO
      ENDDO
      itemp=itemp+nev(itheta)
   ENDDO
ENDDO


DO i=1,nchi*ntheta
   DO j=1,63 !nbig2 tmpmod GregParker
      cc2dtot(i,j,irho)=temparray(i,j)
   ENDDO
ENDDO
ngri=nchi*ntheta


!------------------------------------
!IF needs to be debugged to see the ortho-normalization
!------------------------------------
IF(.false.)THEN
DO i=1,ngri
   DO j=1,ngri
      htemp=0.0d0
      DO k=1,ngri
      htemp=htemp+temparray(k,i)*temparray(k,j)
      ENDDO
      IF(i==j)THEN
           IF(ABS(htemp-1.0d0)>1.0d-6)THEN
           IF(i>ngrj.AND.ABS(htemp)>1.0d-6)THEN
           WRITE(Out_Unit, * )'wrong 2dbasis=',i,j,htemp
           STOP 'error 2D'
           ELSEIF(i<=ngrj)THEN
           WRITE(Out_Unit, * )'wrong 2dbasis=',i,j,htemp
           STOP 'error 2D'
           ENDIF
         ENDIF
      ELSE
           IF(ABS(htemp)>1.0d-6)THEN
         WRITE(Out_Unit, * )'wrong 2dbasis=',i,j,htemp
         STOP 'error 2D'
         ENDIF
      ENDIF
   ENDDO
ENDDO
ENDIF

!--------------------------------------
!WRITE in cc2dtot1 for irho
!--------------------------------------
IF(jtot/=0)THEN
   DO i=1,ntheta*nchi
      DO j=1,i
         WRITE(CC2dtemp_Unit)cc2dtemp(i,j)
      ENDDO
   ENDDO
ENDIF

DEALLOCATE(temparray)
!--------------------------------------

RETURN
END
