SUBROUTINE chctrans(rho, ntheta, nchi, peigv, debug, ltheta, ngood)
USE Numeric_Kinds_Module
USE Numbers_Module
USE FileUnits_Module
USE DVR_Module
USE DVRMasses_Module
USE OneD_Module
USE DVR2_Module
USE Parms_Module
USE TwoD_Module,h2d=>csurf1
!
!     this SUBROUTINE takes the 1-d hamiltonians from hmat1d and the
!     coupling matrices from angmom and turns them into a full 2-d
!     hamiltonian matrix.
!     the resulting matrix is stored in h2d which
!     is THEN diagonalized to give eigenvectors (in h2d) and
!     eigenvalues (in e2d).
!
IMPLICIT NONE
INTEGER ngri, itheta, ntheta, i, j, ichi, ngrj, ntempi, jtheta, ntempj
INTEGER ltheta, i4, item4, i1, item1, i2, nchi, in, k, ngood, ithmode, ibasis
REAL(Kind=WP_Kind) rtwomu, rho, sinth, theta, sinthh, eff, waa, htemp
REAL(Kind=WP_Kind) sd, subd, eigenav
LOGICAL peigv,debug
EXTERNAL mprint
INTEGER, PARAMETER :: ntvmax= nt2max*(nt2max+1)/2
!
!     h2d--2-d hamiltonian
!     later, when h2d is no longer needed, the final 2-d eigenvectors
!     are stored in h2d.
!
!COMMON /twod/ h2d(nt2max,nt2max),       csurf2(ndvptmax,nsfnmax), ovlp(nsfnmax,nsfnmax)
!
COMMON /waste/sd(ntvmax),    subd(nt2max)
rtwomu=One/(twomu*rho*rho)
ngri=0
IF(Ntheta>NthetMax.or.Nchi>NchiMax)THEN
   WRITE(*,*)"Ntheta or Nchi exceeds maximum allowed values"
   WRITE(*,*)"Ntheta=",Ntheta,"  NthetMax=",NthetMax
   WRITE(*,*)"Nchi=",Nchi,"  NchiMax=",NchiMax
   STOP "Chtrans"
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
         cc1d(i,j,itheta)=Zero
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
   sinth=One/(One-pttheta(itheta)**2)
   theta=acos(pttheta(itheta))
   sinthh=One/sin(theta/Two)**2
   !
   !     begin inner loop over theta
   !
   DO jtheta=itheta,ntheta
      ntempj=nev(jtheta)
      eff=Zero
      !
      !     for diagonal elements (both theta points and chi points
      !     are the same), must add in effective potential
      !
      IF(itheta==jtheta)eff=Fifteen/Four-Sixteen*ltheta*ltheta/sin(theta)**2
      waa=rtwomu*(antothet(itheta,jtheta)+eff)
      !
      !     now transform the (diagonal) block of the coupling matrix
      !
      DO i4=1,ntempi
         item4=ngri+i4
         DO i1=1,ntempj
            item1=ngrj+i1
            htemp=Zero
            !
            !     because waa is diagonal in chi, we only need to DO one loop
            !     over chi. this is the heart of the fast transformation
            !
            DO i2=1,nchi
               htemp=htemp+cc1d(i2,i4,itheta)*cc1d(i2,i1,jtheta)*waa
            ENDDO
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
IF(debug) CALL printm(mprint,8,"htot",h2d     ,ngri  ,ngri  ,nt2max)
!
!     diagonalize h2d, the lower triangle of which is stored in the
!     waste array sd. after diagonalization, e2d will contain the
!     eigenvalues and the eigenvectors will be stored columnwise in
!     h2d.
!
k=0
DO i=1,ngri
   DO j=1,i
      k=k+1
      sd(k)=h2d(i,j)
   ENDDO
ENDDO

WRITE(Out_Unit,*)'nt2max=',nt2max,size(h2d),size(h2d,1),size(h2d,2)
DO i=1,nt2max
   DO j=1,nt2max
      h2d(i,j)=Zero
   ENDDO
ENDDO

CALL tdiag_DVR(sd,e2d,h2d,subd,ngri,nt2max,ntvmax)

WRITE(Out_Unit,'(/1x," 2-d h eigenvalues  rho = ",f17.13/)') rho
WRITE(Out_Unit,'(5(1x,f13.9))') (e2d(i),i=1,ngood)

EigenAv=Zero
DO ithmode=1,ngood
   eigenav=eigenav+e2d(ithmode)
ENDDO
eigenav=eigenav/ngood

WRITE(Msg_Unit,*)' rho, eigenav=',rho,eigenav
WRITE(Msg_Unit,*)'(eigenev(ibasis),ibasis=1,MIN(ngood,20)) with ngood=', ngood
WRITE(Msg_Unit,*)(e2d(ibasis),ibasis=1,MIN(ngood,20))
!
IF(peigv) CALL printm(mprint,8,"cc  ",h2d     ,ngri  ,ngri  ,nt2max)
RETURN
ENDSUBROUTINE Chctrans
