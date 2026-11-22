SUBROUTINE Numerov(jval, hstep, mu, v, nmax, eignumerov, Psi, elec, MaxYterms, neig, MaxVibFull)
USE Numeric_Kinds_Module
USE Numbers_Module
USE CriticalDist_Module
USE FileUnits_Module
USE FileUnits_OneDim_Module
USE NArran_Module
USE DiatomicPot_Module
USE OneDim_Parms_Module, Only: jmin, MinElec, IARRAN
USE Convrsns_Module
IMPLICIT NONE
LOGICAL, PARAMETER:: IDbug=.True.
LOGICAL n1, n2
INTEGER nmax, jval, nu, nbound, nquasi, elec, i, numnu, nYterms, info, MaxYterms, neig, NContu, iteration, MaxVibFull
INTEGER, PARAMETER :: lwork=100
REAL(KIND=WP_Kind) work(0:lwork)
REAL(KIND=WP_Kind) hstep, mu, normsquared, Etol, NormInner, NormOuter, Norm
REAL(KIND=WP_Kind) E, Emin, Elower, Eupper, Diff, DiffLower, DiffUpper, DiffMin, EBest, DiffBest, EGuess, EDenom
REAL(KIND=WP_Kind) Psi(0:nmax,0:nmax), v(0:nmax), elist(-1:20000), eignumerov(0:nmax)
REAL(KIND=WP_Kind), ALLOCATABLE :: Yfit(:), amat(:,:)
REAL(KIND=WP_Kind) :: autocm=2.194746313705d5


elist=Vinf(elec)+Two
!CALL TimeInfoSub('Numerov','Begin Numerov')
IF(re.eq.Zero)RETURN
n1=.false.
n2=.false.
Psi=Fifteen
nbound=0
nquasi=0
nu=-1
E=MINVAL(v)    !-de
WRITE(Out_Unit,'(A)')"Numerov Method for Calculating Diatomic Eigenenergies and Wave Functions."
WRITE(Numerov_Unit,*)
WRITE(Numerov_Unit,'(A)')"Numerov Method for Calculating Diatomic Eigenenergies and Wave Functions."
IF(jval.eq.jmin.and.elec.eq.MinElec)THEN
   WRITE(Numerov_Unit,'(A)')"Propagating Solution"
   WRITE(Numerov_Unit,'(A, ES15.7)')"hstep=", hstep
ENDIF
DO WHILE (E<Vinf(elec))
   nu=nu+1
   IF(nu>MaxVibFull)EXIT
   IF(nu==0)THEN
      WRITE(Out_Unit,*)
      WRITE(Out_Unit,*) "Morse Estimate Parameters:"
      WRITE(Out_Unit,*) "De=",De
      WRITE(Out_Unit,*) "Re=",Re
      WRITE(Out_Unit,*) "Beta=",Beta
      WRITE(Out_Unit,*)
      Etol=(EigNumerov(nu+1)-EigNumerov(nu))/Two      !!(EigNumerov(nu)-de)/(Ten**2)
      elower=EigNumerov(nu)-Etol
      eupper=EigNumerov(nu)+Etol
      E=EigNumerov(nu)
   ELSE
      Etol=(EigNumerov(nu)-EigNumerov(nu-1))/Two
      IF(Etol<=Zero)THEN
         nu=nu-1
         EXIT
      ENDIF
      elower=EigNumerov(nu)-Etol
      eupper=EigNumerov(nu)+Etol
      E=EigNumerov(nu)
   ENDIF
   Etol=Eupper-Elower
   CALL NumerOver(nmax, mu, hstep, v, Elower, Psi(:,nu), DiffLower, nu)
   CALL NumerOver(nmax, mu, hstep, v, E,      Psi(:,nu), Diff, nu)
   CALL NumerOver(nmax, mu, hstep, v, Eupper, Psi(:,nu), Diffupper, nu)
   EBest=E
   DiffBest=Diff
   iteration=0
   DO WHILE(ABS(Etol/E).gt.100*Epsilon(One).and.diff.gt.100*Epsilon(One))
      iteration=iteration+1
      !IF(nu==3)THEN
      !   WRITE(numerov_unit,'(A10,I5,6ES15.7)')"iteration=",iteration, Elower, E, Eupper
      !   WRITE(numerov_unit,'(15x,6ES15.7)')DiffLower, DiffBest, DiffUpper
      !   WRITE(numerov_unit,*)
      !ENDIF
      
      IF(MOD(iteration,2)==1)THEN
         Emin=(Elower+E)/Two
      ELSE
         Emin=(Eupper+E)/Two
      ENDIF
      
      CALL NumerOver(nmax, mu, hstep, v, Emin, Psi(:,nu), Diffmin, nu)
      
      IF(Diffmin<DiffBest)THEN
         E=Emin
         EBest=E
         DiffBest=Diffmin
         IF(MOD(iteration,2)==1)THEN
            Elower=(Elower+E)/Two
         ELSE
            Eupper=(Eupper+E)/Two
         ENDIF
         !IF(IDbug.and.nu==3)WRITE(Numerov_Unit,*)"loc3: Elower, DiffLower", Elower, DiffLower
      ELSE
         E=EBest
         EBest=EBest
         DiffBest=DiffBest
         IF(MOD(iteration,2)==1)THEN
            Elower=(Elower+E)/Two
         ELSE
            Eupper=(Eupper+E)/Two
         ENDIF
         !IF(IDbug.and.nu==3)WRITE(Numerov_Unit,*)"loc4: Eupper, DiffUpper", Eupper, DiffUpper
      ENDIF
      Etol=Eupper-Elower
      E=(Eupper+Elower)/Two
      CALL NumerOver(nmax, mu, hstep, v, E, Psi(:,nu), Diff, nu)
   ENDDO
   !IF(ABS(Diff).gt.1.d-6)Exit

   normsquared=hstep*DOT_PRODUCT(Psi(:,nu),Psi(:,nu))
   Psi(0,nmax:nu)=Psi(0,nmax:nu)/SQRT(normsquared)
   IF(nu>=0)elist(nu)=E
   IF(ABS(elist(nu)-elist(nu-1)).lt.1.d-12.or.elist(nu)>=Vinf(elec))EXIT 
   NormInner=hstep*Dot_Product(Psi(0:100,nu),Psi(0:100,nu))
   NormOuter=hstep*Dot_Product(Psi(100:nmax,nu),Psi(100:nmax,nu))
   Norm     =hstep*Dot_Product(Psi(0:nmax,nu),Psi(0:nmax,nu))
   NormInner=NormInner/Norm
   NormOuter=NormOuter/Norm
   IF(E<Vinf(elec))THEN
      IF(E<Vinf(elec))THEN
         IF(nbound==0)THEN
            WRITE(Out_Unit,*)"Ith_Arrangement=",Ith_Arrangement, "  Diatomic Molecule= ",TRIM(DiatomicPot)
            WRITE(Out_Unit,'(A,I0)')"Number of distances=",nmax
            WRITE(Out_Unit,'(a,i3,a,i3)')"Numerov Bound State Energies for jval=", jval," elec=",elec
            WRITE(Numerov_Unit,*)"Ith_Arrangement=",Ith_Arrangement, "  Diatomic Molecule= ",TRIM(DiatomicPot)
            WRITE(Numerov_Unit,'(A,I0)')"Number of distances=",nmax
            WRITE(Numerov_Unit,'(a,i3,a,i3)')"Numerov Bound State Energies for jval=", jval," elec=",elec
         ENDIF
         nbound=nbound+1
         WRITE(Screen_Unit,'(2(a,i4),a,3(es11.3,a),a)')"nu=", nu, " jval=", jval, ' Energy=',e, "(Ha)", e*autoev,"(eV) &
           ",e*autocm,"(cm-1)", ' Bound'
         WRITE(Out_Unit,'(2(a,i4),a,3(es23.15,a),a)')"nu=", nu, " jval=", jval, '  &
		Energy=',E, "(Ha)  ", E*autoev,"(eV)  ",E*autocm,"(cm-1)", ' Bound'
         WRITE(Numerov_Unit,'(2(a,i4),a,3(es23.15,a),a)')"nu=", nu, " jval=", jval, '  &
		Energy=',E, "(Ha)  ", E*autoev,"(eV)  ",E*autocm,"(cm-1)", ' Bound'
      ELSEIF(e.le.v0rc.and.NormInner>NormOuter)THEN
         IF(nquasi==0)THEN
            WRITE(Out_Unit,*)
            WRITE(Out_Unit,'(a,i3,a,i3)')"Numerov Quasibound State Energies for jval=", jval," elec=",elec
            WRITE(Numerov_Unit,*)
            WRITE(Numerov_Unit,'(a,i3,a,i3)')"Numerov Quasibound State Energies for jval=", jval," elec=",elec
         ENDIF
         nquasi=nquasi+1
         WRITE(Screen_Unit,'(2(a,i4),a,3(es11.3,a),a)')"nu=", nu, " jval=", jval, ' Energy=',e, "(Ha)", e*autoev,"(eV) &
           ",e*autocm,"(cm-1)", ' Quasibound'
         WRITE(Out_Unit,'(2(a,i4),a,3(es23.15,a),a)')"nu=", nu, " jval=", jval, &
		'  Energy=',E, "(Ha)  ", E*autoev,"(eV)  ",E*autocm,"(cm-1)", ' Quasibound'
         WRITE(Numerov_Unit,'(2(a,i4),a,3(es23.15,a,a))')"nu=", nu, " jval=", jval, &
		'  Energy=',E, "(Ha)  ", E*autoev,"(eV)  ",E*autocm,"(cm-1)", ' Quasibound'
      ELSE
         IF(NContu.eq.0)THEN
            WRITE(Out_Unit,*)
            WRITE(Out_Unit,'(a,i3,a,i3)')"DAF Continuum State Energies for jval=", jval," elec=",elec
            WRITE(DAF_Unit,*)
            WRITE(DAF_Unit,'(a,i3,a,i3)')"DAF Continuum State Energies for jval=", jval," elec=",elec
         ENDIF
         NContu=NContu+1
         WRITE(Screen_Unit,'(2(a,i4),a,3(es11.3,a),a)')"nu=", nu, " jval=", jval, ' Energy=',e, "(Ha)", e*autoev,"(eV) &
              ",e*autocm,"(cm-1)", ' Continuum'
         WRITE(Out_Unit,'(2(a,i4),a,3(es23.15,a),a)')"nu=", nu, "  jval=", jval, ' Energy=',e, "(Ha)", e*autoev,"(eV) &
              ",e*autocm,"(cm-1)", ' Continuum'
         WRITE(DAF_Unit,'(2(a,i4),a,3(es23.15,a),a)')"nu=", nu, " jval=", jval, ' Energy=',e, "(Ha)", e*autoev,"(eV) &
              ",e*autocm,"(cm-1)", ' Continuum'
      ENDIF
   ENDIF
   IF(E<Vinf(elec))THEN
      WRITE(NumerovBin_Unit)nu, jval, elec
      WRITE(NumerovBin_Unit)elist(nu)
      WRITE(NumerovBin_Unit)Psi(0:nmax,nu)
   ENDIF
ENDDO
IF(nquasi>0)THEN
   WRITE(Numerov_Unit,*)"Error, nquasi:", nquasi
   !stop "nquasi>0"
ENDIF
numnu=nbound+nquasi
neig=nbound
WRITE(NumerovBin_Unit)-numnu, jval, elec

nYterms=min(MaxYterms,numnu)
ALLOCATE(Yfit(numnu))
ALLOCATE(Amat(numnu,nYterms+1))
DO nu=0,numnu-1
   Yfit(nu+1)=elist(nu)
   DO i=0,nYterms-1
      amat(nu+1,i+1)=(nu+Half)**i
   ENDDO
ENDDO
IF(Numnu.gt.0)CALL dgels('N',numnu,nYterms,1,amat,numnu,Yfit,numnu,work,lwork,info)
WRITE(Out_Unit,*)
WRITE(Out_Unit,*) "Numerov Results"
WRITE(Out_Unit,*)'Info=',Info
WRITE(Out_Unit,*) "Spectroscopic Constants: Y(nu) for j=", jval, "and elec=",elec
WRITE(Spectra_Unit,*)
WRITE(Spectra_Unit,*) "Numerov Results"
WRITE(Spectra_Unit,'(a,i3,a,i3,a,i3)')"elec=",elec,"  jval=",jval,"  nYterms=",nYterms
DO i=0,nYterms-1
   WRITE(Out_Unit,'(a,i2,a,ES16.7,"(Ha)",ES16.7,"(cm^-1)")')"Y(",i,")=",Yfit(i+1),Yfit(i+1)/cmm1toau
   WRITE(Spectra_Unit,'(a,i3,a,ES23.14,"(Ha)",ES23.14,"(cm^-1)")')"  i=",i,"  Yfit=",Yfit(i+1),Yfit(i+1)/cmm1toau
ENDDO
IF(jval==0)THEN
   !WRITE(Out_Unit,'(2A6,5A17,",",5A17)')"elec,","jval,","Yfit(0)(Ha),","Yfit(1)(Ha),","Yfit(2)(Ha),","Yfit(3)(Ha),","Yfit(4)(Ha),","Yfit(0)(cm^-1),","Yfit(1)(cm^-1),","Yfit(2)(cm^-1),","Yfit(3)(cm^-1),","Yfit(4)(cm^-1),"
   !WRITE(Spectra_Unit,'(2A6,5A17,",",5A17)')"elec,","jval,","Yfit(0)(Ha),","Yfit(1)(Ha),","Yfit(2)(Ha),","Yfit(3)(Ha),","Yfit(4)(Ha),","Yfit(0)(cm^-1),","Yfit(1)(cm^-1),","Yfit(2)(cm^-1),","Yfit(3)(cm^-1),","Yfit(4)(cm^-1),"
   WRITE(YfitPltNumerov_Unit,'(2A6,5A17,",",5A17)')"elec,","jval,","Yfit(0)(Ha),","Yfit(1)(Ha),","Yfit(2)(Ha),","Yfit(3)(Ha),","Yfit(4)(Ha),","Yfit(0)(cm^-1),","Yfit(1)(cm^-1),","Yfit(2)(cm^-1),","Yfit(3)(cm^-1),","Yfit(4)(cm^-1),"
ENDIF
WRITE(YfitPltNumerov_Unit,'(2(I5,","),5(ES16.7,","),",",5(ES16.7,","))')elec, jval,(Yfit(i+1),i=0,nYterms-1),(Yfit(i+1)/cmm1toau,i=0,nYterms-1)
WRITE(YfitBin_Unit)elec,jval,numnu,nbound,nquasi,nYTerms
WRITE(YfitBin_Unit)"Numerov    ",DiatomicPot, Iarran, elec, jval, (Yfit(i+1),i=0,nYterms-1)

WRITE(NumerovBin_Unit)elec,jval,numnu,nbound,nquasi,nYterms
WRITE(NumerovBin_Unit)(yfit(i),i=1,nYterms)
DEALLOCATE(Yfit,AMAT)

!CALL TimeInfoSub('Numerov','  End Numerov')

!STOP "TmpMODGregParker Numerov"

RETURN
END SUBROUTINE Numerov
