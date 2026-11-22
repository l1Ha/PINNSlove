SUBROUTINE FourierGrid(n, Pot, mu, eigenvalues, Psi, Hamil, hstep, neig, MaxVibFull)
   USE Numeric_Kinds_Module
   USE Numbers_Module
   USE CriticalDist_Module
   USE FileUnits_Module
   USE FileUnits_OneDim_Module
   USE QstateDAF_Module
   USE OneDim_Parms_Module
   USE NArran_Module
   USE DiatomicPot_Module
   USE Convrsns_Module
   IMPLICIT NONE  
   CHARACTER(LEN=1), PARAMETER:: Jobz='V', Range='V', UpLo='U'
   ! PORT FIX (2026-10-01): lwork was only assigned on the first call
   ! (jval==jmin .and. elec==MinElec) and was not SAVEd, so later calls passed
   ! an undefined LWORK to DSYEV ("parameter number 8 had an illegal value").
   INTEGER, SAVE :: lwork = 1000000
   INTEGER info, NBound, NQuasi, NContu, Numnu, Neig, iwork(1000000), ifail(1000000)
   INTEGER i, j, n, l, npoints, nyterms, nup, MaxVibFull
   REAL(KIND=wp_kind) prefac, hstep, mu, e, sum, pot0save, NormInner, NormOuter, Norm, Const1, Const2
   REAL(KIND=wp_kind) Pot(0:2*n), eigenvalues(0:2*n), Psi(0:2*n,0:2*n), Hamil(0:2*n,0:2*n)
   REAL(KIND=wp_kind) s_work(1000000)
   REAL(KIND=WP_Kind), ALLOCATABLE:: work(:), q(:,:), amat(:,:), yFit(:), DVR(:,:)
   REAL(KIND=WP_Kind) :: autocm=2.194746313705d5
   
   
  ! CALL TimeInfoSub('Numerov','Begin FourierGrid')
   WRITE(Out_Unit,'(A)')"FourierGrid Method for Calculating Diatomic Eigenenergies and Wave Functions."
   WRITE(FourierGrid_Unit,*)
   WRITE(FourierGrid_Unit,'(A)')"FourierGrid Method for Calculating Diatomic Eigenenergies and Wave Functions."
   pot0save=Pot(0)
   IF(ABS(Pot(0))>ABS(Pot(1)))Pot(0)=Pot(1)

!     F. Gogtas, G.G. Balint-Kurti and C.C. Marston, QCPE Program No. 647, (1993).
!     The Fourier Grid Hamiltonian method is used in the program.
!     This method is fully described in:
!     C.C. Marston and G.G.Balint-Kurti, J.Chem. Phys.,91,3571(1989).
!     and
!     G.G. Balint-Kurti, R.N. Dixon and C.C. Marston, Internat. Rev. Phys. Chem.,11, 317 (1992).
!
!     The program uses an even number of grid points. The Hamiltonian
!     matrix elements are calculated using analytic formula described
!     in the second of the above references.
!     The analytical Hamiltonian expression given in the reference
!     contains a small error.  The formula should read: 
!
!                H(i,j) = {(h**2)/(4*m*(L**2)} *
!                          [(N-1)(N-2)/6 + (N/2)] + V(Xi),   if i=j
!
!                H(i,j) = {[(-1)**(i-j)] / m } *
!                          { h/[2*L*sin(pi*(i-j)/N)]}**2 ,     if i#j
!
!     The eigenvalues of the Hamiltonian matrix which lie below the
!     asymptotic (large x) value of V(x) correspond to the bound state
!     energies. The corresponding eigenvectors are the representation
!     of the bound state wavefunctions on a regular one dimensional grid.
Hamil=Zero
npoints=2*n+1
Const1=(Pi**2)/(mu*((npoints*hstep)**2))
Const2=((Npoints-1)*(Npoints-2)/Six + Npoints/Two)
DO i=0,npoints-1
   DO j=0,npoints-1
      IF(i==j)THEN
         Hamil(i,j) = Const1*Const2 + Pot(i)
      ELSE
         Hamil(i,j) = ((-1)**(i-j))*Const1/(SIN(Pi*(i-j)/Npoints))**2
      ENDIF
   ENDDO
ENDDO
     
IF(jval.eq.jmin.and.elec.eq.MinElec)THEN
   lwork=1000000
   WRITE(FourierGrid_Unit,*)"Calling dsyev"
   WRITE(FourierGrid_Unit,*)"lwork=",lwork
ENDIF
CALL dsyev(Jobz, UpLo, npoints, Hamil, npoints, eigenvalues, s_work, lwork, info)
neig=npoints
!CALL dsyevx(Jobz,Range, UpLo, npoints, Hamil, npoints, v0re, vasy, 1, 1, 1.d-12, neig, eigenvalues, Psi, npoints, s_work, lwork, iwork, ifail, info)
IF(info/=0)THEN
   WRITE(FourierGrid_Unit,*)"Error returned from dsyevx: Info=", Info
   STOP "Error in FourierGrid when dsyevx was called"
ENDIF
   
NBound=0
NQuasi=0
NContu=0
neig=npoints
neig=MaxVibFull
Psi=Hamil
DO nu=0,Neig
   Psi(0,nu)=Zero
   Psi(npoints-1,nu)=Zero
   e=eigenvalues(nu)
   NormInner=hstep*Dot_Product(Psi(0:100,nu),Psi(0:100,nu))
   NormOuter=hstep*Dot_Product(Psi(100:npoints-1,nu),Psi(100:npoints-1,nu))
   Norm     =hstep*Dot_Product(Psi(0:npoints-1,nu),Psi(0:npoints-1,nu))
   NormInner=NormInner/Norm
   NormOuter=NormOuter/Norm
   IF(e<vasy)THEN
      IF(NBound.eq.0)THEN
         WRITE(Out_Unit,*)"Ith_Arrangement=",Ith_Arrangement, "  Diatomic Molecule= ",TRIM(DiatomicPot)
         WRITE(Out_Unit,'(A,I0)')"Number of distances=",npoints
         WRITE(Out_Unit,'(a,i3,a,i3)')"FourierGrid Bound State Energies for jval=", jval," elec=",elec
         WRITE(FourierGrid_Unit,*)"Ith_Arrangement=",Ith_Arrangement, "  Diatomic Molecule= ",TRIM(DiatomicPot)
         WRITE(FourierGrid_Unit,'(A,I0)')"Number of distances=",npoints
         WRITE(FourierGrid_Unit,'(a,i3,a,i3)')"FourierGrid Bound State Energies for jval=", jval," elec=",elec
      ENDIF
      NBound=NBound+1
      WRITE(Screen_Unit,'(2(a,i4),a,3(es11.3,a),a)')"nu=", nu, " jval=", jval, ' Energy=',e, "(Ha)", e*autoev,"(eV) &
           ",e*autocm,"(cm-1)", ' Bound'
      WRITE(Out_Unit,'(2(a,i4),a,3(es23.15,a),a)')"nu=", nu, " jval=", jval, ' Energy=',e, "(Ha)", e*autoev,"(eV) &
           ",e*autocm,"(cm-1)", ' Bound'
      WRITE(FourierGrid_Unit,'(2(a,i4),a,3(es23.15,a),a)')"nu=", nu, " jval=", jval, ' Energy=',e, "(Ha)", e*autoev,"(eV) &
           ",e*autocm,"(cm-1)",   ' Bound'
   ELSEIF(e.le.v0rc.and.NormInner>NormOuter)THEN
      IF(NQuasi.eq.0)THEN
         WRITE(Out_Unit,*)
         WRITE(Out_Unit,'(a,i3,a,i3)')"FourierGrid Quasibound State Energies for jval=", jval," elec=",elec
         WRITE(FourierGrid_Unit,*)
         WRITE(FourierGrid_Unit,'(a,i3,a,i3)')"FourierGrid Quasibound State Energies for jval=", jval," elec=",elec
      ENDIF
      NQuasi=NQuasi+1
      WRITE(Screen_Unit,'(2(a,i4),a,3(es11.3,a),a)')"nu=", nu, " jval=", jval, ' Energy=',e, "(Ha)", e*autoev,"(eV) &
           ",e*autocm,"(cm-1)", ' Quasibound'
      WRITE(Out_Unit,'(2(a,i4),a,3(es23.15,a),a)')"nu=", nu, "  jval=", jval, ' Energy=',e, "(Ha)", e*autoev,"(eV) &
           ",e*autocm,"(cm-1)", ' Quasibound'
      WRITE(FourierGrid_Unit,'(2(a,i4),a,3(es23.15,a),a)')"nu=", nu, " jval=", jval, ' Energy=',e, "(Ha)", e*autoev,"(eV) &
           ",e*autocm,"(cm-1)", ' Quasibound'
   ELSE
      IF(NContu.eq.0.and.continuum)THEN
         WRITE(Out_Unit,*)
         WRITE(Out_Unit,'(a,i3,a,i3)')"FourierGrid Continuum State Energies for jval=", jval," elec=",elec
         WRITE(FourierGrid_Unit,*)
         WRITE(FourierGrid_Unit,'(a,i3,a,i3)')"FourierGrid Continuum State Energies for jval=", jval," elec=",elec
      ENDIF
      IF(continuum)THEN
         NContu=NContu+1
         WRITE(Screen_Unit,'(2(a,i4),a,3(es11.3,a),a)')"nu=", nu, " jval=", jval, ' Energy=',e, "(Ha)", e*autoev,"(eV) &
           ",e*autocm,"(cm-1)", ' Continuum'
         WRITE(Out_Unit,'(2(a,i4),a,3(es23.15,a),a)')"nu=", nu, "  jval=", jval, ' Energy=',e, "(Ha)", e*autoev,"(eV) &
           ",e*autocm,"(cm-1)", ' Continuum'
         WRITE(FourierGrid_Unit,'(2(a,i4),a,3(es23.15,a),a)')"nu=", nu, " jval=", jval, ' Energy=',e, "(Ha)", e*autoev,"(eV) &
           ",e*autocm,"(cm-1)", ' Continuum'
      ENDIF
   ENDIF
ENDDO

Psi=Psi/SQRT(hstep)
numnu=nbound+nquasi
nYterms=min(MaxYterms,numnu)
ALLOCATE(Yfit(numnu))
ALLOCATE(Amat(numnu,nYterms))
WRITE(FourierGridBin_Unit)elec,jval,numnu,nbound,nquasi,nYTerms
DO nu=0,numnu-1
   Yfit(nu+1)=eigenvalues(nu)
   WRITE(FourierGridBin_Unit)nu,eigenvalues(nu)
   WRITE(FourierGridBin_Unit)psi(:,nu)
   DO i=0,nYterms-1
      amat(nu+1,i+1)=(nu+Half)**i
   ENDDO
ENDDO

IF(Numnu.ne.0)THEN
   ALLOCATE(DVR(0:numnu-1,0:numnu-1))
   DO nu=0,numnu-1
      DO nup=0,numnu-1
         SUM=Zero
         DO I=0,2*n
            SUM=SUM+psi(i,nu)*psi(i,nup)
         ENDDO
         DVR(nu,nup)=SUM
      ENDDO
   ENDDO
   DVR=DVR*hstep
   DEALLOCATE(DVR)
ENDIF


IF(.not.ALLOCATED(work))ALLOCATE(work(lwork))
IF(Numnu.gt.0)CALL dgels('N',numnu,nYterms,1,amat,numnu,Yfit,numnu,work,lwork,info)
WRITE(Out_Unit,*)
WRITE(Out_Unit,*) "FourierGrid Results"
WRITE(Out_Unit,*)'Info=',Info
WRITE(Out_Unit,*) "Spectroscopic Constants: Y(nu) for j=", jval, "and elec=",elec
WRITE(Spectra_Unit,*)
WRITE(Spectra_Unit,*) "FourierGrid Results"
WRITE(Spectra_Unit,'(a,i3,a,i3,a,i3)')"elec=",elec,"  jval=",jval,"  nYterms=",nYterms
DO i=0,nYterms-1
   WRITE(Out_Unit,'(a,i2,a,ES16.7,"(Ha)",ES16.7,"(cm^-1)")')"Y(",i,")=",Yfit(i+1),Yfit(i+1)/cmm1toau
   WRITE(Spectra_Unit,'(a,i3,a,ES23.14,"(Ha)",ES23.14,"(cm^-1)")')"  i=",i,"  Yfit=",Yfit(i+1),Yfit(i+1)/cmm1toau
ENDDO
IF(jval==0)THEN
   !WRITE(Out_Unit,'(2A6,5A17,",",5A17)')"elec,","jval,","Yfit(0)(Ha),","Yfit(1)(Ha),","Yfit(2)(Ha),","Yfit(3)(Ha),","Yfit(4)(Ha),","Yfit(0)(cm^-1),","Yfit(1)(cm^-1),","Yfit(2)(cm^-1),","Yfit(3)(cm^-1),","Yfit(4)(cm^-1),"
   !WRITE(Spectra_Unit,'(2A6,5A17,",",5A17)')"elec,","jval,","Yfit(0)(Ha),","Yfit(1)(Ha),","Yfit(2)(Ha),","Yfit(3)(Ha),","Yfit(4)(Ha),","Yfit(0)(cm^-1),","Yfit(1)(cm^-1),","Yfit(2)(cm^-1),","Yfit(3)(cm^-1),","Yfit(4)(cm^-1),"
   WRITE(YfitPltFourier_Unit,'(2A6,5A17,",",5A17)')"elec,","jval,","Yfit(0)(Ha),","Yfit(1)(Ha),","Yfit(2)(Ha),","Yfit(3)(Ha),","Yfit(4)(Ha),","Yfit(0)(cm^-1),","Yfit(1)(cm^-1),","Yfit(2)(cm^-1),","Yfit(3)(cm^-1),","Yfit(4)(cm^-1),"
ENDIF
WRITE(YfitPltFourier_Unit,'(2(I5,","),5(ES16.7,","),",",5(ES16.7,","))')elec, jval,(Yfit(i+1),i=0,nYterms-1),(Yfit(i+1)/cmm1toau,i=0,nYterms-1)
WRITE(YfitBin_Unit)elec,jval,numnu,nbound,nquasi,nYTerms
WRITE(YfitBin_Unit)"FourierGrid",DiatomicPot, Iarran, elec, jval, (Yfit(i+1),i=0,nYterms-1)

WRITE(FourierGridBin_Unit)(Yfit(i),i=1,nYterms)
DEALLOCATE(Yfit,AMAT)

neig=nbound+nquasi
!CALL TimeInfoSub('FourierGrid','  End FourierGrid')

   
 Pot(0)=pot0save
RETURN
ENDSUBROUTINE FourierGrid