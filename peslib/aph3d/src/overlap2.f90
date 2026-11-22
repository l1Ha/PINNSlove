SUBROUTINE overlap2(rholst, rhonow, smat, check, nchan, ifirst)
USE Numeric_Kinds_Module
USE FileUnits_Module
USE pbasis_Module
IMPLICIT NONE
LOGICAL little, medium, full, ifirst
INTEGER nchan, ithcll, ithsub, lam, nmodes, j, iinit, jj, ifinal, i
INTEGER megamax, n, np
REAL(Kind=WP_Kind) rho1, rho2, rholst, rhonow, olast, xlarge
REAL(Kind=WP_Kind)  smat(nchan,nchan), check(nchan)
DATA ithcll/0/,ithsub/0/
DATA little/.false./,medium/.false./,full/.false./

!revised by X.L.
megamax=maxmega
IF(Little)WRITE(Out_Unit,*)'In Overlap2'
CALL popt('overlap2', little, medium, full, ithcll, ithsub)
!-----------------------------------------------------------------------
! Skip over the previous lambda values.
!-----------------------------------------------------------------------
DO lam = megamin, lammin-1
   READ(ovr_unit) nmodes,rho1,rho2
   IF (abs(rholst-rho1)/rho1>1.d-7)THEN
      WRITE(Out_Unit,110) rholst,rho1
      STOP 'overlap2'
   ENDIF
   DO j=1,nmodes
         READ(ovr_unit)
   ENDDO
   READ(ovr_unit)
ENDDO

iinit=1
jj=1
ifinal=0
DO 5 lam = lammin, lammax
!          IF (iregion/='aph    ')THEN
!-----------------------------------------------------------------------
! store the old eigenfunctions in chiold
!-----------------------------------------------------------------------
!             nsum=nnuj(1)+nnuj(2)+nnuj(3)
!             DO 10 i=1,nsum
!    10       chiold(i)=chinuj(i)
!-----------------------------------------------------------------------
! calculate the new eigenfunctions
!-----------------------------------------------------------------------
!             CALL upsiln (vecold,tstore,xsq,xk,rhonow,chinuj)
!             DO 20 np=1,nchan
!             DO 20 n=1,nchan
!    20       smat(n,np)=0.0d0
!-----------------------------------------------------------------------
! calculate the overlaps.
!-----------------------------------------------------------------------
!temp notice this requires that all of the nhermt's be exactly the
!temp same for all arrangement channels.
!DO 50 i=1,nhermt(1)
!DO 30 ichan=1,3
!IF (c(ichan)==0.0d0.or..NOT.integrat(ichan)) go to 30
!cnow=c(ichan)/rhonow
!bnow=b(ichan)/rhonow
!clast=c(ichan)/rholst
!blast=b(ichan)/rholst
!IF(scheme==1)THEN
!   ynow=xpth(i,ichan)
!   ylast=(cnow*xpth(i,ichan)+bnow-blast)/clast
!   CALL dhep (h(1,ichan),ylast,noscil)
!   expfac(ichan)=exp((ynow**2-ylast**2)/2)
!ENDIF
!30 CONTINUE
!DO 50 n=1,nchan
!ichan=kchan(n)
!jl=jrot3(n)
!nvib=mvib3(n)
!kin=nvib*noscil(ichan)
!cold=0.0d0
!IF(scheme==2)THEN
!   karran=ichan
!   ylast=(cnow*xpth(i,ichan)+bnow-blast)/clast
!   CALL hwave (h(1,ichan),ylast,jl,noscil(ichan))
!   expfac(ichan)=1.0d0
!ENDIF
!DO k=1,noscil(ichan)
!   cold=cold+chlast(kin+k,jl,ichan)*h(k,ichan)
!ENDDO
!cold=cold*expfac(ichan)
!DO 50 np=1,nchan
!nujmp=(nuj3(np)-1)*nhermt(ichan)
!IF (kchan(n)/=kchan(np)) go to 50
!IF (lorb3(n)/=lorb3(np)) go to 50
!IF (jrot3(n)/=jrot3(np)) go to 50
!inp=i+nujmp
!smat(n,np)=smat(n,np)+cold*chinuj(inp)* wpth(i,ichan)*sqrt(rholst/rhonow)
!50 CONTINUE
!ELSE
      IF(lam/=lammin) READ(ovr_unit)
      READ(ovr_unit,err=180) nmodes,rho1,rho2
      IF(little)WRITE(Out_Unit,100) nmodes,rho1,rho2
      IF (abs(rholst-rho1)/rho1>1.d-7)THEN
         WRITE(Out_Unit,110) rholst,rho1
         STOP 'overlap2'
      ENDIF
      IF (nmodes<nlam(lam))THEN
         WRITE(Out_Unit,130) nmodes,nlam(lam)
         STOP 'overlap2'
      ENDIF
      ifinal=ifinal+nlam(lam)
      DO j=1,nmodes
         IF(j<=nlam(lam))THEN
            READ(ovr_unit,err=170) (smat(i,jj),i=iinit,ifinal)
            jj=jj+1
         ELSE
            READ(ovr_unit)
         ENDIF
      ENDDO
      iinit=iinit+nlam(lam)

!           ENDIF
!-----------------------------------------------------------------------
! print the overlap2 matrix
!-----------------------------------------------------------------------
IF (little) WRITE(Out_Unit,140) rholst,rhonow
IF (full) CALL MxOut(smat,nchan,nchan)
!
IF (little.AND.ifirst)THEN
   !-----------------------------------------------------------------------
   ! find the largest of the off-diagonal elements.
   !-----------------------------------------------------------------------
   xlarge=0.0d0
   DO n=2,nchan
      DO np=1,n-1
         IF (xlarge<abs(smat(n,np))) xlarge=abs(smat(n,np))
      ENDDO
   ENDDO
   IF(little)WRITE(Out_Unit,150) xlarge

   !-----------------------------------------------------------------------
   ! calculate olast.
   !-----------------------------------------------------------------------
   olast=1.0d0

   DO n=1,nchan
      IF (olast>abs(smat(n,n))) olast=abs(smat(n,n))
      IF (smat(n,n)==0.0d0) go to 80
      IF (olast>1.0d0/abs(smat(n,n))) olast=1.0d0/abs(smat(n,n))
   80 CONTINUE
   ENDDO

   IF (little) WRITE(Out_Unit,160) olast
   CALL basisck (smat, check, nchan)
   olast=1.0d0

   DO n=1,nchan
      IF (olast>abs(check(n))) olast=abs(check(n))
      IF (check(n)==0.0d0) go to 81
      IF (olast>1.0d0/abs(check(n))) olast=1.0d0/abs(check(n))
   81 CONTINUE
   ENDDO
   WRITE(Out_Unit,160) olast
ENDIF


!DO i=1,nchan
!   xsq(i)=xksq3(i)
!   xk(i)=sqrt(abs(xsq(i)))
!ENDDO
!           ENDIF
5     CONTINUE

!-----------------------------------------------------------------------
! Skip over the remaining lambda values.
!-----------------------------------------------------------------------
DO lam = lammax+1, megamax
   READ(ovr_unit)
   READ(ovr_unit) nmodes,rho1,rho2
   IF (abs(rholst-rho1)/rho1>1.d-7)THEN
      WRITE(Out_Unit,110) rholst,rho1
      STOP 'overlap2'
   ENDIF
   DO j=1,nmodes
         READ(ovr_unit)
   ENDDO
ENDDO 
IF(Little)WRITE(Out_Unit,*)'Leaving Oerlap2'
RETURN
!
100 FORMAT (//,'nmodes=',i5,'  rho1=',f10.5,'  rho2=',f10.5)
110 FORMAT (//,'***error in Overlap2***: rholst/=rho1',2f10.5)
130 FORMAT (//,'***error in Overlap2*** nmodes<nlam(lam)',2i5)
140 FORMAT (//,' overlap2 matrix for rho-last=',f10.5,'  rhonow=',f10.5)
150 FORMAT (//,' the magnitude of the largest off-diagonal',' overlap2 matrix element is ',f10.5)
160 FORMAT (//,'magnitude of the smallest diagonal overlap2',f10.5)
!
170   WRITE(msg_unit,*)' Error: error reading S_Matrix from Ovrlp_ABM in overlap2'
WRITE(msg_unit,*)' iinit, ifinal=', iinit, ifinal
STOP 'overlap2'
180   WRITE(msg_unit,*)' Error: error reading nmodes,rho1,rho2 from',' Ovrlp_ABM'
STOP 'overlap2'
ENDSUBROUTINE Overlap2
