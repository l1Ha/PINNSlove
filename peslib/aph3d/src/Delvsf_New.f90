SUBROUTINE Delvsf_New (nt, sinu, cosu, dsinu, dcosu, edsinu, edcosu, dsinued, dcosued, fk, oopsln)
!
!----------------------------------------------------------------------
! calculates the space-fixed delves asymptotic solutions.
! written by g. a. parker
! ----------------------------------------------------------------------
USE Numbers_Module
USE FileUnits_Module
USE nstate_Module
USE Narran_Module
USE parms_Module
USE region_Module
USE VFunc_Module
USE NumNuj_Module
USE thetas_Module
USE Arrch_Module
USE Qall_Module
USE gaussb_MODULE
USE GaussQuady_Module
USE Integrat_Module
USE Oops_Module
USE Rhonum_Module
USE Spectro_Module
USE opts_Module
USE Masses_Module
USE etachn_Module
USE das_Module
USE AzBzCz_Module
USE kmat_MODULE
IMPLICIT NONE
INTEGER NT, IthCll, IthSub, Nsum, Ichan, k, jacob, nvib, jl, lorb, kin
INTEGER idelv, kdelve
REAL(KIND=WP_Kind) kmthke, kmhke, oopstemp, cg, bg, cf, bf, yf
REAL(KIND=WP_Kind) ThetaD, STheta, CTheta, STrans, yg, gamf, fjacob
REAL(KIND=WP_Kind) sfj, a, zj, dzj, zy, dzy, chi, z, dchi, aj, ay, cj, cy
REAL(KIND=WP_Kind) dj, dy, edaj, eday, edcj, edcy, eddj, eddy
LOGICAL little,medium,full
REAL(Kind=WP_Kind) fk(nt), oopsln(1), h(100,3), expfac(3)
REAL(Kind=WP_Kind) sinu(nt,nt), cosu(nt,nt), dsinu(nt,nt), dcosu(nt,nt)
!----------------------------------------------------------------------
REAL(Kind=WP_Kind) edsinu(nt,nt), edcosu(nt,nt), dsinued(nt,nt), dcosued(nt,nt)
REAL(KIND=WP_Kind) thetd_temp(maxhermt,narran), bfaci, az1, bz1, cz1
REAL(KIND=WP_Kind) xpthtemp(maxhermt)
!----------------------------------------------------------------------
! Determine printing options.
!----------------------------------------------------------------------
DATA ithcll/0/,ithsub/0/
DATA little/.false./,medium/.false./,full/.false./
CALL popt('Delvsf_New  ', little, medium, full, ithcll, ithsub)

newkmat=.true.
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
! PORT (2026-10-11): long-range matching.  upsiln resets rhoj=rhoval (the
! deck's $rhovalue, i.e. the Delves-region end at 8.1 a0).  When the
! R-matrix has been propagated out to APH3D_RF by Jacobi_LongRange_Prop,
! every geometric factor of the asymptotic matching must be evaluated at
! that SAME hyperradius: strans = rhoj*cos(theta_d) (the Jacobi
! separation), the Delves theta grid (basis_delves at rho=rhonum), the
! Morse coordinate map yg, the Jacobian bfaci and the Bessel argument
! z = k*rhoj*cos(theta_d).  A hybrid convention (Bessel argument at one
! radius, geometry at another) corrupts the per-arrangement phase
! structure of the extracted S-matrix.
BLOCK
   CHARACTER(LEN=32) :: env_rf_m
   INTEGER :: stat_rf_m
   REAL(Kind=WP_Kind) :: rf_m
   CALL get_environment_variable('APH3D_RF', env_rf_m, status=stat_rf_m)
   IF (stat_rf_m == 0) THEN
      READ(env_rf_m, *) rf_m
      IF (rf_m > rhoj) THEN
         rhoj   = rf_m
         rhonum = rf_m
      END IF
   END IF
END BLOCK
!-----------------------------------------------------------------------
! initialize all matrices to zero.
!-----------------------------------------------------------------------
CALL vsets (nt*nt, sinu, 1, Zero)
CALL vsets (nt*nt, cosu, 1, Zero)
CALL vsets (nt*nt, dsinu, 1, Zero)
CALL vsets (nt*nt, dcosu, 1, Zero)
CALL vsets (nt*nt, edsinu, 1, Zero)
CALL vsets (nt*nt, edcosu, 1, Zero)
CALL vsets (nt*nt, dsinued, 1, Zero)
CALL vsets (nt*nt, dcosued, 1, Zero)
!-----------------------------------------------------------------------
! determine the asymptotic delves solutions
!-----------------------------------------------------------------------
DO ichan = 1,narran
   IF((c(ichan)/=Zero).AND.integrat(ichan))THEN
   cg=c(ichan)
   bg=b(ichan)
   cf=c(ichan)/rhoj
   bf=b(ichan)/rhoj
   karran=ichan
   CALL basis_delves(xpth(1,ichan),thetd_temp,az,bz,cz,ichan)
   az1=az(ichan)
   bz1=bz(ichan)
   cz1=cz(ichan)
   CALL basis_jacobi(xpth(1,ichan),az(ichan),bz(ichan),cz(ichan))
   DO k=1,nhermt(ichan)
      yf=xpth(k,ichan)
! start my change we need the new method to determine everything.
! We do need to start with the theta_delves
! then we need to know the S, and s
! revised by X.L. 01-23-2007
! step 1: theta_delves
! new way is
      thetad=thetd_temp(k,ichan)
      stheta=sin(thetad)
      ctheta=cos(thetad)
      strans=rhoj*ctheta
      yg=az(ichan)*(rhoj*sin(thetad))-bz(ichan)/(rhoj*sin(thetad))+cz(ichan)
      !WRITE(Out_Unit,*)k,ichan,stheta,thetad
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
             STOP 'Delvsf_New'
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
             STOP 'Delvsf_New'
         ENDIF
!!!               expfac(ichan) = exp(gamf*(yf-yg)/2)*sqrt(rhoj*cf/cg)
         expfac(ichan) = 1.d0
      ELSE
         WRITE(Out_Unit,*)'Error: scheme must be 1 or 2'
         WRITE(Out_Unit,*)'scheme=',scheme
         STOP 'Delvsf_New'
      ENDIF
      DO jacob=1,nt
         IF(kchan(jacob)==ichan)THEN
         nvib=mvib3(jacob)
         jl=jrot3(jacob)
         fjacob=fk(jacob)
         lorb=lorb3(jacob)
         sfj=sqrt(fjacob)
         ! PORT (2026-10-11): the long-range matching radius is now carried
         ! consistently by rhoj (set from APH3D_RF in Kmatrx_Main), so the
         ! Bessel argument uses the proper theta-dependent Jacobi separation
         ! z = k * rhoj * cos(theta_d) — identical to the original convention.
         z=fjacob*strans
         IF(xksq3(jacob)>0.0)THEN
            CALL rbes (lorb,z,zj,dzj,zy,dzy)
         ELSE
            CALL kbes (lorb,z,zj,dzj,zy,dzy)
         ENDIF
         kin=nvib*(3*noscil(ichan))
!-----------------------------------------------------------------------
! make the jacobi vibrational function, chi, and it's first
! derivative, dchi, at the delves point yg
!-----------------------------------------------------------------------
         karran=ichan
         CALL JaVib_New(yg, 3*noscil(ichan), h, jl, expfac, kin, chi, dchi,ichan, cg,rhoj*stheta,rhoj)
         DO idelv=1,nt
            IF((jrot3(idelv)==jrot3(jacob)).AND.(lorb3(idelv)==lorb3(jacob)).AND.(kchan(idelv)==kchan(jacob)))THEN
!-----------------------------------------------------------------------
! project the asymptotic jacobi solutions onto the delves basis
!-----------------------------------------------------------------------
            kdelve=(nuj3(idelv)-1)*nhermt(ichan)+k
            oopstemp=oopsln(kdelve)
            bfaci=sqrt(rhoj)/(az1/ctheta/ctheta+bz1/stheta/stheta)
            aj=wpth(k,ichan)*oopstemp*chi*zj/sfj*bfaci
            ay=wpth(k,ichan)*oopstemp*chi*zy/sfj*bfaci
            cj=wpth(k,ichan)*oopstemp*chi*ctheta*dzj*sfj*bfaci
            cy=wpth(k,ichan)*oopstemp*chi*ctheta*dzy*sfj*bfaci
            dj=wpth(k,ichan)*oopstemp*stheta*dchi*zj/sfj*bfaci
            dy=wpth(k,ichan)*oopstemp*stheta*dchi*zy/sfj*bfaci
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
            edaj = wpth(k,ichan)*oopstemp*chi*kmthke*(z*dzj - 0.5d0*zj)*bfaci
            eday = wpth(k,ichan)*oopstemp*chi*kmthke*(z*dzy - 0.5d0*zy)*bfaci
            edcj = wpth(k,ichan)*oopstemp*chi*ctheta*kmhke*(0.5d0*dzj + (z)**(-1)*(lorb*(lorb + 1) - z**2)*zj)*bfaci
            edcy = wpth(k,ichan)*oopstemp*chi*ctheta*kmhke*(0.5d0*dzy + (z)**(-1)*(lorb*(lorb + 1) - z**2)*zy)*bfaci
            eddj = wpth(k,ichan)*oopstemp*stheta*dchi*kmthke*(z*dzj - 0.5d0*zj)*bfaci
            eddy = wpth(k,ichan)*oopstemp*stheta*dchi*kmthke*(z*dzy - 0.5d0*zy)*bfaci
!-----------------------------------------------------------------------
!   determine the energy derivatives of sinu, cosu, dsinu, and dcosu
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
! PORT (2026-10-11): fingerprint of the matching geometry and the projected
! free-wave matrices, so the actual matching convention of a run can be
! verified from the log alone.
BLOCK
   REAL(Kind=WP_Kind) :: fp_e, fp_f
   fp_e = SUM(ABS(sinu(1:nt,1:nt)))
   fp_f = SUM(ABS(cosu(1:nt,1:nt)))
   WRITE(Out_Unit,'(A,ES23.15,A,ES23.15)') &
      ' DBG Delvsf_New: rhoj=', rhoj, '  fpE=', fp_e
   WRITE(Out_Unit,'(A,ES23.15)') ' DBG Delvsf_New: fpF=', fp_f
   WRITE(*,'(A,ES23.15,A,ES23.15)') &
      ' DBG Delvsf_New: rhoj=', rhoj, '  fpE=', fp_e
END BLOCK
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
ENDSUBROUTINE Delvsf_New
