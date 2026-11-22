SUBROUTINE Delvsf_Old (nt, sinu,   cosu,   dsinu,   dcosu, edsinu, edcosu, dsinued, dcosued, fk, oopsln)
!
!----------------------------------------------------------------------
! calculates the space-fixed delves asymptotic solutions.
! written by g. a. parker
! ----------------------------------------------------------------------
USE FileUnits_Module
USE nstate_Module
USE Narran_Module
USE parms_Module
USE region_Module
USE VFunc_Module
USE NumNuj_Module
USE thetas_Module
USE Oops_Module
USE Qall_Module
USE GaussQuady_Module
USE Integrat_Module
USE Arrch_Module
USE Spectro_Module
USE opts_Module
USE Masses_Module
USE etachn_Module
USE das_Module
IMPLICIT NONE
INTEGER ithcll, ithsub, nsum, ichan, k, nvib, jl, lorb, idelv, kdelve, jacob
INTEGER kin, nt
REAL(Kind=WP_Kind) kmthke, kmhke, oopstemp, cg, bg, cf, bf, yf, thetad, stheta
REAL(Kind=WP_Kind) ctheta, strans, yg, gamf, fjacob, sfj, z, zj, dzj, chi
REAL(Kind=WP_Kind) zy, dzy, dchi, aj, ay, cj, dy, edaj, cy, dj, eday
REAL(Kind=WP_Kind) edcy, eddy, edcj, eddj
LOGICAL little,medium,full
REAL(Kind=WP_Kind) fk(nt), oopsln(1), h(100,3), expfac(3)
REAL(Kind=WP_Kind) sinu(nt,nt), cosu(nt,nt), dsinu(nt,nt), dcosu(nt,nt)
!----------------------------------------------------------------------
REAL(Kind=WP_Kind) edsinu(nt,nt), edcosu(nt,nt), dsinued(nt,nt), dcosued(nt,nt)
!----------------------------------------------------------------------
! Determine printing options.
!----------------------------------------------------------------------
DATA ithcll/0/,ithsub/0/
DATA little/.false./,medium/.false./,full/.false./
CALL popt('Delvsf_Old  ', little, medium, full, ithcll, ithsub)

! ----------------------------------------------------------------------
! Calculate the delves eigenfuctions.
! ----------------------------------------------------------------------
ioops = .true.
iregion = 'devles '
CALL upsiln (sinu, cosu, dsinu, fk, rhoj, chinuj)
! ----------------------------------------------------------------------
! store the delves eigenfunctions in oopsln
! ----------------------------------------------------------------------
nsum=nnuj(1)+nnuj(2)+nnuj(3)
CALL vcopy (nsum, oopsln, chinuj)
!-----------------------------------------------------------------------
! calculate the jacobi eigenfunctions and store them in chinuj
!-----------------------------------------------------------------------
ioops=.false.
iregion='jacobi '
CALL upsiln (sinu, cosu, dsinu, fk, rhoj, chinuj)
!-----------------------------------------------------------------------
! initialize all matrices to zero.
!-----------------------------------------------------------------------
CALL vsets (nt*nt, sinu, 1, 0.0d0)
CALL vsets (nt*nt, cosu, 1, 0.0d0)
CALL vsets (nt*nt, dsinu, 1, 0.0d0)
CALL vsets (nt*nt, dcosu, 1, 0.0d0)
CALL vsets (nt*nt, edsinu, 1, 0.0d0)
CALL vsets (nt*nt, edcosu, 1, 0.0d0)
CALL vsets (nt*nt, dsinued, 1, 0.0d0)
CALL vsets (nt*nt, dcosued, 1, 0.0d0)
!-----------------------------------------------------------------------
! determine the asymptotic delves solutions
!-----------------------------------------------------------------------
DO ichan = 1,narran
   IF((c(ichan)/=0.0).AND.integrat(ichan))THEN
   cg=c(ichan)
   bg=b(ichan)
   cf=c(ichan)/rhoj
   bf=b(ichan)/rhoj
   DO k=1,nhermt(ichan)
      yf=xpth(k,ichan)
      thetad=cf*yf+bf
      stheta=sin(thetad)
      ctheta=cos(thetad)
      strans=rhoj*ctheta
      yg=(rhoj*stheta-bg)/cg
      IF(scheme==1)THEN
         IF((yf**2-yg**2)/2 > 700.d0)THEN
             WRITE(Out_Unit,*)'Error: Overflow will occur'
             WRITE(Out_Unit,*)'(yf**2-yg**2)/2=',(yf**2-yg**2)/2
             WRITE(Out_Unit,*)'yf,yg=',yf,yg
             WRITE(Out_Unit,*)'k,ichan',k,ichan
             WRITE(Out_Unit,*)'cg,bg,cf,bf=',cg,bg,cf,bf
             WRITE(Out_Unit,*)'rhoj=',rhoj
             WRITE(Out_Unit,*)'thetad=',thetad
             WRITE(Out_Unit,*)'stheta,ctheta,strans=',stheta,ctheta,strans
             STOP 'Delvsf_Old'
         ENDIF
         expfac(ichan)=exp((yf**2-yg**2)/2)*sqrt(rhoj*cf/cg)
      ELSEIF(scheme==2)THEN
         gamf = 2.d0*eta(ichan)*usys/dscale(ichan)
         IF(gamf*(yf-yg)/2 > 700.d0)THEN
             WRITE(Out_Unit,*)'Error: Overflow will occur'
             WRITE(Out_Unit,*)'(yf-yg)/2=',(yf-yg)/2
             WRITE(Out_Unit,*)'yf,yg=',yf,yg
             WRITE(Out_Unit,*)'k,ichan',k,ichan
             WRITE(Out_Unit,*)'cg,bg,cf,bf=',cg,bg,cf,bf
             WRITE(Out_Unit,*)'rhoj=',rhoj
             WRITE(Out_Unit,*)'thetad=',thetad
             WRITE(Out_Unit,*)'stheta,ctheta,strans=',stheta,ctheta,strans
             WRITE(Out_Unit,*)'gamf=',gamf
             STOP 'Delvsf_Old'
         ENDIF
!!!               expfac(ichan) = exp(gamf*(yf-yg)/2)*sqrt(rhoj*cf/cg)
         expfac(ichan) = 1.d0
      ELSE
         WRITE(Out_Unit,*)'Error: scheme must be 1 or 2'
         WRITE(Out_Unit,*)'scheme=',scheme
         STOP 'Delvsf_Old'
      ENDIF
      DO jacob=1,nt
         IF(kchan(jacob)==ichan)THEN
         nvib=mvib3(jacob)
         jl=jrot3(jacob)
         fjacob=fk(jacob)
         lorb=lorb3(jacob)
         sfj=sqrt(fjacob)
         z=fjacob*strans
         IF(xksq3(jacob)>0.0)THEN
            CALL rbes (lorb,z,zj,dzj,zy,dzy)
         ELSE
            CALL kbes (lorb,z,zj,dzj,zy,dzy)
         ENDIF
         kin=nvib*noscil(ichan)
!-----------------------------------------------------------------------
! make the jacobi vibrational function, chi, and it's first
! derivative, dchi, at the delves point yg
!-----------------------------------------------------------------------
         karran=ichan
         CALL JaVib_Old(yg, noscil(ichan), h, jl, expfac, kin, chi, dchi, ichan, cg)
         DO idelv=1,nt
            IF((jrot3(idelv)==jrot3(jacob)).AND.(lorb3(idelv)==lorb3(jacob)).AND.(kchan(idelv)==kchan(jacob)))THEN
!-----------------------------------------------------------------------
! project the asymptotic jacobi solutions ono the delves basis
!-----------------------------------------------------------------------
            kdelve=(nuj3(idelv)-1)*nhermt(ichan)+k
            oopstemp=oopsln(kdelve)
            aj=wpth(k,ichan)*oopstemp*chi*zj/sfj
            ay=wpth(k,ichan)*oopstemp*chi*zy/sfj
            cj=wpth(k,ichan)*oopstemp*chi*ctheta*dzj*sfj
            cy=wpth(k,ichan)*oopstemp*chi*ctheta*dzy*sfj
            dj=wpth(k,ichan)*oopstemp*stheta*dchi*zj/sfj
            dy=wpth(k,ichan)*oopstemp*stheta*dchi*zy/sfj
            sinu(idelv,jacob)=sinu(idelv,jacob)+aj
            cosu(idelv,jacob)=cosu(idelv,jacob)+ay
            dsinu(idelv,jacob)=dsinu(idelv,jacob)+(aj/rhoj/2+cj+dj)
            dcosu(idelv,jacob)=dcosu(idelv,jacob)+(ay/rhoj/2+cy+dy)
!-----------------------------------------------------------------------
!   obtain kmthke = (k**-3/2)*kE (e-deriv. of k)
!                 = (fjacob**-5/2)*usys
!   obtain kphke  = k*kpmhke
!-----------------------------------------------------------------------
            kmthke = fjacob**(-5.d0/2.d0)*usys
            kmhke  = fjacob*kmthke
            kmhke  = fjacob**(-3.d0/2.d0)*usys
!-----------------------------------------------------------------------
!   determine the energy derivatives of aj, yj, cj, cy, dj, and dy
!-----------------------------------------------------------------------
            edaj = wpth(k,ichan)*oopstemp*chi*kmthke*(z*dzj - 0.5d0*zj)
            eday = wpth(k,ichan)*oopstemp*chi*kmthke*(z*dzy - 0.5d0*zy)
            edcj = wpth(k,ichan)*oopstemp*chi*ctheta*kmhke*(0.5d0*dzj + (z)**(-1)*(lorb*(lorb + 1) - z**2)*zj)
            edcy = wpth(k,ichan)*oopstemp*chi*ctheta*kmhke*(0.5d0*dzy + (z)**(-1)*(lorb*(lorb + 1) - z**2)*zy)
            eddj = wpth(k,ichan)*oopstemp*stheta*dchi*kmthke*(z*dzj - 0.5d0*zj)
            eddy = wpth(k,ichan)*oopstemp*stheta*dchi*kmthke*(z*dzy - 0.5d0*zy)
!-----------------------------------------------------------------------
!                determine the energy derivatives of sinu, cosu, dsinu, and dcosu
!-----------------------------------------------------------------------
            edsinu(idelv,jacob) = edsinu(idelv,jacob) + edaj
            edcosu(idelv,jacob) = edcosu(idelv,jacob) + eday
            dsinued(idelv,jacob) = dsinued(idelv,jacob) + edaj/rhoj/2 + edcj + eddj
            dcosued(idelv,jacob) = dcosued(idelv,jacob) + eday/rhoj/2 + edcy + eddy
            ENDIF
         ENDDO
         ENDIF
      ENDDO
   ENDDO
   ENDIF
ENDDO
!-----------------------------------------------------------------------
! WRITE out the projected matrics.
!-----------------------------------------------------------------------
IF(full)THEN
    WRITE(Out_Unit,70)
    CALL MxOutL(sinu, nt, nt, 0, 'space', 'space')
    WRITE(Out_Unit,80)
    CALL MxOutL(cosu, nt, nt, 0, 'space', 'space')
    WRITE(Out_Unit,90)
    CALL MxOutL(dsinu, nt, nt, 0, 'space', 'space')
    WRITE(Out_Unit,100)
    CALL MxOutL(dcosu, nt, nt, 0, 'space', 'space')
ELSEIF(medium)THEN
    WRITE(Out_Unit,70)
    CALL MxOut(sinu, nt, nt)
    WRITE(Out_Unit,*)'Energy Derivative'
    CALL MxOut(edsinu, nt, nt)
    WRITE(Out_Unit,80)
    CALL MxOut(cosu, nt, nt)
    WRITE(Out_Unit,*)'Energy Derivative'
    CALL MxOut(edcosu, nt, nt)
    WRITE(Out_Unit,90)
    CALL MxOut(dsinu, nt, nt)
    WRITE(Out_Unit,*)'Energy Derivative'
    CALL MxOut(dsinued, nt, nt)
    WRITE(Out_Unit,100)
    CALL MxOut(dcosu, nt, nt)
    WRITE(Out_Unit,*)'Energy Derivative'
    CALL MxOut(dcosued, nt, nt)
ENDIF

RETURN
!
70 FORMAT(//,'sinu-matrix')
80 FORMAT(//,'cosu-matrix')
90 FORMAT(//,'dsinu-matrix')
100 FORMAT(//,'dcosu-matrix')
ENDSUBROUTINE Delvsf_Old
