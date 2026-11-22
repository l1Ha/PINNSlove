SUBROUTINE overlap(rholst,rhonow,smat,nchan,vecold,tstore,xsq,xk,chiold,olast)
USE Numeric_Kinds_Module
USE FileUnits_MODULE
USE parms_MODULE
USE nstate_MODULE
USE Narran_Module
USE thetas_MODULE
USE region_MODULE
USE Qall_Module
USE GaussQuady_Module
USE NumNuj_Module
USE VFunc_Module
USE cherm1_MODULE
USE cherm2_MODULE
USE schem_MODULE
USE Integrat_Module
USE Arrch_Module
IMPLICIT NONE
SAVE
INTEGER ithcll, ithsub, nsum, i, np, n, ichan, jl, nvib, k, nujmp, nmodes
INTEGER kin, nchan, inp, j
LOGICAL little,medium,full
REAL(Kind=WP_Kind) rhonow, cnow, bnow, clast, rholst, blast, ynow, ylast, cold
REAL(Kind=WP_Kind) rho1, rho2, xlarge, olast
REAL(Kind=WP_Kind) chiold(*)   ! PORT: assumed-size (sequence association)
REAL(Kind=WP_Kind) smat(nchan,nchan), tstore(*), vecold(*)   ! PORT: assumed-size
REAL(Kind=WP_Kind) xsq(*), xk(*)   ! PORT: assumed-size
REAL(Kind=WP_Kind) h(100,3), hdn(100,3)
REAL(Kind=WP_Kind) expfac(3), check(nstate)
DATA ithcll/0/,ithsub/0/
DATA little/.false./,medium/.false./,full/.false./
CALL popt('overlap ', little, medium, full, ithcll, ithsub)
little=.TRUE.
IF(iregion/='aph    ')THEN
   !-----------------------------------------------------------------------
   ! store the old eigenfunctions in chiold
   !-----------------------------------------------------------------------
   nsum=nnuj(1)+nnuj(2)+nnuj(3)
   DO i=1,nsum
      chiold(i)=chinuj(i)
   ENDDO
   !-----------------------------------------------------------------------
   ! calculate the new eigenfunctions
   !-----------------------------------------------------------------------
   CALL upsiln (vecold,tstore,xsq,xk,rhonow,chinuj)
   DO np=1,nchan
      DO n=1,nchan
         smat(n,np)=0.0d0
      ENDDO
   ENDDO

   !-----------------------------------------------------------------------
   ! calculate the overlaps.
   !-----------------------------------------------------------------------
   !temp notice this requires that all of the nhermt's be exactly the
   !temp same for all arrangement channels.
   DO i=1,nhermt(1)
      ! PORT FIX (2026-10-02): this loop precomputes PER-ARRANGEMENT quantities
      ! (c, b, integrat, xpth, noscil, h, expfac are all dimensioned narran),
      ! so it must run over the 3 arrangements, not over the nchan channels
      ! (channels are mapped to arrangements below via kchan(n)).
      ! Previously it read past those arrays (nchan=40 > narran=3) -> crash.
      DO ichan=1,narran
         IF(c(ichan)/=0.0d0.AND.integrat(ichan))THEN
            cnow=c(ichan)/rhonow
            bnow=b(ichan)/rhonow
            clast=c(ichan)/rholst
            blast=b(ichan)/rholst
            IF(scheme==1)THEN
               ynow=xpth(i,ichan)
               ylast=(cnow*xpth(i,ichan)+bnow-blast)/clast
               CALL dhep (h(1,ichan),ylast,3*noscil(ichan))
               expfac(ichan)=exp((ynow**2-ylast**2)/2)
            ENDIF
         ENDIF
      ENDDO
      DO n=1,nchan
         ichan=kchan(n)
         jl=jrot3(n)
         nvib=mvib3(n)
         kin=nvib*3*noscil(ichan)
         cold=0.0d0
         IF(scheme==2)THEN
            karran=ichan
            ylast=(cnow*xpth(i,ichan)+bnow-blast)/clast
            CALL hwave (h(1,ichan),ylast,jl,3*noscil(ichan))
            expfac(ichan)=1.0d0
         ENDIF
         DO k=1,3*noscil(ichan)
            cold=cold+chlast(kin+k,jl,ichan)*h(k,ichan)
         ENDDO
         cold=cold*expfac(ichan)
         DO np=1,nchan
            nujmp=(nuj3(np)-1)*nhermt(ichan)
            IF((kchan(n)==kchan(np)).AND.(lorb3(n)==lorb3(np)).AND.(jrot3(n)==jrot3(np)))THEN
               inp=i+nujmp
               smat(n,np)=smat(n,np)+cold*chinuj(inp)*wpth(i,ichan)*sqrt(rholst/rhonow)
            ENDIF
         ENDDO
      ENDDO
   ENDDO
ELSE

   !WRITE(Out_Unit,*)'reading from Ovr_Unit pos1'
   READ(Ovr_Unit) nmodes,rho1,rho2

   !IF(little)WRITE(Out_Unit,100) nmodes,rho1,rho2
   IF(ABS(rholst-rho1)/rho1>1.d-7)THEN
      WRITE(Out_Unit,110) rholst,rho1
   ENDIF

   IF(nmodes<nchan)THEN
      WRITE(Out_Unit,130) nmodes,nchan
      WRITE(Msg_Unit,130) nmodes,nchan
      STOP 'Stopping in overlap'
   ENDIF

   !WRITE(Out_Unit,*)'reading from Ovr_Unit pos2'
   DO j=1,nmodes
      IF(j<=nchan) READ(Ovr_Unit) (smat(i,j),i=1,nchan)
      IF(j>nchan) READ(Ovr_Unit)
   ENDDO
   !WRITE(Out_Unit,*)'end reading Ovr_Unit pos2'
ENDIF

!-----------------------------------------------------------------------
! print the overlap matrix
!-----------------------------------------------------------------------
IF(little) WRITE(Out_Unit,140) rholst,rhonow
IF(full) CALL MxOut(smat,nchan,nchan)
!-----------------------------------------------------------------------
! find the largest of the off-diagonal elements.
!-----------------------------------------------------------------------
xlarge=0.0d0
DO n=2,nchan
   DO np=1,n-1
      IF(xlarge<ABS(smat(n,np))) xlarge=ABS(smat(n,np))
   ENDDO
ENDDO
IF(little)WRITE(Out_Unit,150) xlarge

!-----------------------------------------------------------------------
! calculate olast.
!-----------------------------------------------------------------------
olast=1.0d0
DO n=1,nchan
   IF(olast>ABS(smat(n,n))) olast=ABS(smat(n,n))
   IF(smat(n,n)/=0.0d0)THEN
      IF(olast>1.0d0/ABS(smat(n,n))) olast=1.0d0/ABS(smat(n,n))
   ENDIF
ENDDO

IF(little) WRITE(Out_Unit,160) olast
CALL basisck (smat, check, nchan)

olast=1.0d0
DO n=1,nchan
   IF(olast>ABS(check(n))) olast=ABS(check(n))
   IF(check(n)/=0.0d0)THEN
      IF(olast>1.0d0/ABS(check(n))) olast=1.0d0/ABS(check(n))
   ENDIF
ENDDO
IF(little) WRITE(Out_Unit,160) olast

DO i=1,nchan
   xsq(i)=xksq3(i)
   xk(i)=sqrt(ABS(xsq(i)))
ENDDO
WRITE(Out_Unit,*)'Returning from overlap'
WRITE(Out_Unit,*)

RETURN
!
100 FORMAT(1x,'nmodes=',i5,'  rho1=',f10.5,'  rho2=',f10.5)
110 FORMAT(1x,'***error in overlap***: rholst/=rho1',2f10.5)
130 FORMAT(1x,'***error in overlap*** nmodes<nlam(lam)',2i5)
140 FORMAT(1x,'overlap matrix for rho-last=',f10.5,'  rhonow=',f10.5)
150 FORMAT(1x,'the magnitude of the largest off-diagonal',' overlap matrix element is ',f10.5)
160 FORMAT(1x,'magnitude of the smallest diagonal overlap',f10.5)
ENDSUBROUTINE OVERLAP
