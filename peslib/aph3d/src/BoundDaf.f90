SUBROUTINE BoundDAF(hstep, mu, v, eig, psi, kinetic, toeplitz, dafx, ab, neig, MaxVibFull)
USE Numeric_Kinds_Module
USE Numbers_Module
USE CriticalDist_Module
USE FileUnits_Module
USE FileUnits_OneDim_Module
USE OneDim_Parms_Module
USE QstateDAF_Module
USE NArran_Module
USE DiatomicPot_Module
USE Convrsns_Module
IMPLICIT NONE
SAVE
CHARACTER(LEN=8), PARAMETER, DIMENSION(7):: Methods(7)=(/'dsyev   ', 'dsyevd  ', 'dsbev   ', 'dsbevd  ', 'dsyevx  ', 'dsyevr  ', 'dsbevx  '/)
CHARACTER(LEN=6), PARAMETER:: Method='dsbevx'
CHARACTER(LEN=1), PARAMETER:: Jobz='V', UpLo='L', Range='V'
CHARACTER(LEN=5) StateType(0:nmax)
CHARACTER(LEN=15) :: Identifier
INTEGER, PARAMETER:: ldab=42,nxval=3
INTEGER kmax, NBound, NQuasi, i, nup, info, j, Neig, lda, lwork, liwork, ldq, NContu, kmaxsave, MaxVibFull
! PORT FIX (2026-10-01): IFAIL of DSYEVX/DSBEVX is an array of length N; it was
! declared as a scalar, so LAPACK wrote N integers over adjacent memory
! (SIGSEGV inside DSBEVX).  nmax is a module constant by call time.
INTEGER ifail(0:nmax)
INTEGER Numnu, nYTerms, s_iwork(1)
INTEGER, ALLOCATABLE:: iwork(:), isuppz(:)
LOGICAL, PARAMETER:: idbug=.false.
Logical Banded
REAL(KIND=WP_Kind), PARAMETER:: AbsTol=1.d-12
REAL(KIND=WP_Kind) mu, hstep, e, sum, NormInner, NormOuter, Norm
REAL(KIND=WP_Kind) v(0:nmax), xval(nxval)
REAL(KIND=WP_Kind) kinetic(0:nmax,0:nmax), toeplitz(0:nmax), dafx(0:nmax,0:ndaf), hermit(0:mval+ndaf+1), ab(0:ldab,0:nmax+ldab)
REAL(KIND=WP_Kind) psi(0:nmax,0:nmax), eig(0:nmax), s_work(100000)
REAL(KIND=WP_Kind), ALLOCATABLE:: work(:), q(:,:), amat(:,:), yFit(:), DVR(:,:)
REAL(KIND=WP_Kind) :: autocm=2.194746313705d5

CALL TimeInfoSub('BoundDAF','Begin BoundDAF')
WRITE(Out_Unit,'(A)')"DAF Method for Calculating Diatomic Eigenenergies and Wave Functions."
WRITE(DAF_Unit,*)
WRITE(DAF_Unit,'(A)')"DAF Method for Calculating Diatomic Eigenenergies and Wave Functions."
kinetic=Zero
ab=Zero
ifail=0
IF(jval.eq.jmin.and.elec.eq.MinElec)THEN
   Identifier="_BoundDaf_"
    CALL get_daf(nxval, xval, hstep, mval, dafx, hermit, ndaf, nmax, toeplitz, Two*mu, kmax, Identifier)
    kmaxsave=kmax
ENDIF
kmax=kmaxsave

IF((Method.eq.Methods(3).or.Method.eq.Methods(4).or.Method.eq.Methods(7)).and.jval.eq.jmin)THEN
   Banded=.true.
   lda=MAX(1,SIZE(ab,1))
   !WRITE(DAF_Unit,*)'LDA=',lda
ELSE IF(jval.eq.jmin)THEN
   Banded=.false.
   lda=MAX(1,SIZE(Kinetic,1))
   !WRITE(DAF_Unit,*)'LDA=',lda
ENDIF
    
CALL kin_radial(nmax+1, toeplitz, kinetic, kmax)

ab=Zero
DO i=0,nmax
   kinetic(i,i)=kinetic(i,i)+v(i)
   DO j=0,i
      IF(j.le.i.and.i.le.min(nmax,j+ldab))THEN
         ab(i-j,j)=kinetic(i,j)
      ENDIF
   ENDDO
ENDDO

IF(idbug)THEN
   WRITE(Dbug_Unit,*)'Hamiltonian Matrix'
   CALL mxout(kinetic, nmax+1, nmax+1)
ENDIF

IF(Method.eq.Methods(1))THEN
   IF(jval.eq.jmin.and.elec.eq.MinElec)THEN
      WRITE(DAF_Unit,*)"Calling dsyev"
      lwork=-1
      CALL dsyev(Jobz, UpLo, nmax+1, kinetic, lda, eig, s_work, lwork, info)
      lwork=s_work(1)
      WRITE(DAF_Unit,*)"lwork=",lwork
      psi=kinetic
      IF(.not.ALLOCATED(work))ALLOCATE(work(lwork))
   ENDIF
   CALL dsyev(Jobz, UpLo, nmax+1, kinetic, lda, eig, work, lwork, info)
   Neig=nmax+1
ELSE IF(Method.eq.Methods(2))THEN
   IF(jval.eq.jmin.and.elec.eq.MinElec)THEN
      WRITE(DAF_Unit,*)"Calling dsyevd"
      lwork=-1
      liwork=-1
      CALL dsyevd(Jobz, UpLo, nmax+1, kinetic, lda, eig, s_work, lwork, s_iwork, liwork, info)
      lwork=s_work(1)
      liwork=s_iwork(1)
      psi=kinetic
      WRITE(DAF_Unit,*)"lwork=",lwork," liwork=",liwork
      IF(.not.ALLOCATED(work))ALLOCATE(work(lwork),iwork(liwork))
   ENDIF
   CALL dsyevd(Jobz, UpLo, nmax+1, kinetic, lda, eig, work, lwork, iwork, liwork, info)
   Neig=nmax+1
ELSE IF(Method.eq.Methods(3))THEN
   IF(jval.eq.jmin.and.elec.eq.MinElec)THEN
      WRITE(DAF_Unit,*)"Calling dsbev"
      lwork=MAX(1,3*(nmax+1)-2)
      WRITE(DAF_Unit,*)"lwork=",lwork
      ALLOCATE(work(lwork))
   ENDIF
   CALL dsbev(Jobz, UpLo, nmax+1, kmax, ab, lda, eig, psi, nmax+1, work, info)
   Neig=nmax+1
ELSE IF(Method.eq.Methods(4))THEN
   IF(jval.eq.jmin.and.elec.eq.MinElec)THEN
      WRITE(DAF_Unit,*)"Calling dsbevd"
      lwork=-1
      liwork=-1
      CALL dsbevd(Jobz, UpLo, nmax+1, kmax, ab, lda, eig, psi, nmax+1, s_work, lwork, s_iwork, liwork, info)
      lwork=s_work(1)
      liwork=s_iwork(1)
      WRITE(DAF_Unit,*)"lwork=",lwork," liwork=",liwork
      ALLOCATE(work(lwork),iwork(liwork))
   ENDIF
   CALL dsbevd(Jobz, UpLo, nmax+1, kmax, ab, lda, eig, psi, nmax+1, work, lwork, iwork, liwork, info)
   Neig=nmax+1
ELSE IF(Method.eq.Methods(5))THEN
   IF(jval.eq.jmin.and.elec.eq.MinElec)THEN
      WRITE(DAF_Unit,*)"Calling dsyevx"
      lwork=-1
      liwork=MAX(1,5*(nmax+1))
      CALL dsyevx(Jobz, Range, UpLo, nmax+1, kinetic, nmax+1, v0re, v0rc, 1, 1, AbsTol, Neig, eig, &
		psi, nmax+1, s_work, lwork, iwork, ifail, info)
      lwork=s_work(1)
      WRITE(DAF_Unit,*)"lwork=",lwork," liwork=",liwork
      ALLOCATE(work(lwork),iwork(liwork))
   ENDIF
   CALL dsyevx(Jobz, Range, UpLo, nmax+1, kinetic, nmax+1, v0re, v0rc, 1, 1, AbsTol, Neig, eig, & 
	psi, nmax+1, work, lwork, iwork, ifail, info)
ELSE IF(Method.eq.Methods(6))THEN
   IF(jval.eq.jmin.and.elec.eq.MinElec)THEN
      WRITE(DAF_Unit,*)"Calling dsyevr"
      ALLOCATE(isuppz(2*MAX(1,nmax+1)))
      lwork=-1
      liwork=-1
      CALL dsyevr(Jobz, Range, UpLo, nmax+1, kinetic, nmax+1, v0re, v0rc, 1, 1, AbsTol, Neig, &
	eig, psi, nmax+1, isuppz, s_work, lwork, s_iwork, liwork, info)
      lwork=s_work(1)
      liwork=s_iwork(1)
      WRITE(DAF_Unit,*)"lwork=",lwork," liwork=",liwork
      ALLOCATE(work(lwork),iwork(liwork))
   ENDIF
   CALL dsyevr(Jobz, Range, UpLo, nmax+1, kinetic, nmax+1, v0re, v0rc, 1, 1, AbsTol, Neig, &
	eig, psi, nmax+1, isuppz, work, lwork, iwork, liwork, info)
ELSE IF(Method.eq.Methods(7))THEN
   IF(jval.eq.jmin.and.elec.eq.MinElec)THEN
      WRITE(DAF_Unit,*)"Calling dsbevx"
      lwork=MAX(1,7*(nmax+1))
      liwork=MAX(1,5*(nmax+1))
      ldq=nmax+1
      WRITE(DAF_Unit,*)"lwork=",lwork," liwork=",liwork," ldq=",ldq
      IF(.not.ALLOCATED(work))ALLOCATE(work(lwork),iwork(liwork),q(0:nmax,0:nmax))
   ENDIF
   CALL dsbevx(Jobz, Range, UpLo, nmax+1, kmax, ab, lda, q, ldq, v0re, v0rc, 1, 1, &
  	AbsTol, Neig, eig, psi, nmax+1, work, iwork, ifail, info)
ELSE
   WRITE(Out_unit,*)'Error Method=', Method
   WRITE(Msg_unit,*)'Error Method=', Method
   WRITE(DAF_unit,*)'Error Method=', Method
   STOP 'BoundDaf'
ENDIF

IF(ANY(ifail/=0).or.info.ne.0)THEN
   WRITE(Msg_Unit,*)"ifail=",ifail,"  info=",info
   WRITE(DAF_Unit,*)"ifail=",ifail,"  info=",info
   STOP "BoundDAF"
ENDIF

NBound=0
NQuasi=0
NContu=0
neig=MaxVibFull
DO nu=0,Neig
   e=eig(nu)
   NormInner=hstep*Dot_Product(Psi(0:100,nu),Psi(0:100,nu))
   NormOuter=hstep*Dot_Product(Psi(100:nmax,nu),Psi(100:nmax,nu))
   Norm     =hstep*Dot_Product(Psi(0:nmax,nu),Psi(0:nmax,nu))
   NormInner=NormInner/Norm
   NormOuter=NormOuter/Norm
   IF(e<vasy)THEN
      IF(NBound.eq.0)THEN
         WRITE(Out_Unit,*)"Ith_Arrangement=",Ith_Arrangement, "  Diatomic Molecule= ",TRIM(DiatomicPot)
         WRITE(Out_Unit,'(A,I0)')"Number of distances=",nmax
         WRITE(Out_Unit,'(a,i3,a,i3)')"DAF Bound State Energies for jval=", jval," elec=",elec
         WRITE(DAF_Unit,*)"Ith_Arrangement=",Ith_Arrangement, "  Diatomic Molecule= ",TRIM(DiatomicPot)
         WRITE(DAF_Unit,'(A,I0)')"Number of distances=",nmax
         WRITE(DAF_Unit,'(a,i3,a,i3)')"DAF Bound State Energies for jval=", jval," elec=",elec
      ENDIF
      NBound=NBound+1
      WRITE(Screen_Unit,'(2(a,i4),a,3(es11.3,a),a)')"nu=", nu, " jval=", jval, ' Energy=',e, "(Ha)", e*autoev,"(eV) &
           ",e*autocm,"(cm-1)", ' Bound'
      WRITE(Out_Unit,'(2(a,i4),a,3(es23.15,a),a)')"nu=", nu, " jval=", jval, ' Energy=',e, "(Ha)", e*autoev,"(eV) &
           ",e*autocm,"(cm-1)", ' Bound'
      WRITE(DAF_Unit,'(2(a,i4),a,3(es23.15,a),a)')"nu=", nu, " jval=", jval, ' Energy=',e, "(Ha)", e*autoev,"(eV) &
           ",e*autocm,"(cm-1)",   ' Bound'
   ELSEIF(e.le.v0rc.and.NormInner>NormOuter)THEN
      IF(NQuasi.eq.0)THEN
         WRITE(Out_Unit,*)
         WRITE(Out_Unit,'(a,i3,a,i3)')"DAF Quasibound State Energies for jval=", jval," elec=",elec
         WRITE(DAF_Unit,*)
         WRITE(DAF_Unit,'(a,i3,a,i3)')"DAF Quasibound State Energies for jval=", jval," elec=",elec
      ENDIF
      NQuasi=NQuasi+1
      WRITE(Screen_Unit,'(2(a,i4),a,3(es11.3,a),a)')"nu=", nu, " jval=", jval, ' Energy=',e, "(Ha)", e*autoev,"(eV) &
           ",e*autocm,"(cm-1)", ' Quasibound'
      WRITE(Out_Unit,'(2(a,i4),a,3(es23.15,a),a)')"nu=", nu, "  jval=", jval, ' Energy=',e, "(Ha)", e*autoev,"(eV) &
           ",e*autocm,"(cm-1)", ' Quasibound'
      WRITE(DAF_Unit,'(2(a,i4),a,3(es23.15,a),a)')"nu=", nu, " jval=", jval, ' Energy=',e, "(Ha)", e*autoev,"(eV) &
           ",e*autocm,"(cm-1)", ' Quasibound'
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
ENDDO

Psi=Psi/SQRT(hstep)
numnu=nbound+nquasi
nYterms=min(MaxYterms,numnu)
ALLOCATE(Yfit(numnu))
ALLOCATE(Amat(numnu,nYterms))
elec=0  !TmpModGregParker
WRITE(DAFBin_Unit)elec,jval,numnu,nbound,nquasi,nYTerms
DO nu=0,numnu-1
   Yfit(nu+1)=eig(nu)
   WRITE(DAFBin_Unit)nu,eig(nu)
   WRITE(DAFBin_Unit)psi(:,nu)
   DO i=0,nYterms-1
      amat(nu+1,i+1)=(nu+Half)**i
   ENDDO
ENDDO

IF(Numnu.ne.0)THEN
   ALLOCATE(DVR(0:numnu-1,0:numnu-1))
   DO nu=0,numnu-1
      DO nup=0,numnu-1
         SUM=Zero
         DO I=0,nmax
            SUM=SUM+psi(i,nu)*psi(i,nup)
         ENDDO
         DVR(nu,nup)=SUM
      ENDDO
   ENDDO
   DVR=DVR*hstep
   DEALLOCATE(DVR)
ENDIF

IF(Numnu.gt.0)CALL dgels('N',numnu,nYterms,1,amat,numnu,Yfit,numnu,work,lwork,info)
WRITE(Out_Unit,*)
WRITE(Out_Unit,*) "DAF Results"
WRITE(Out_Unit,*)'Info=',Info
WRITE(Out_Unit,*) "Spectroscopic Constants: Y(nu) for j=", jval, "and elec=",elec
WRITE(Spectra_Unit,*)
WRITE(Spectra_Unit,*) "DAF Results"
WRITE(Spectra_Unit,'(a,i3,a,i3,a,i3)')"elec=",elec,"  jval=",jval,"  nYterms=",nYterms
DO i=0,nYterms-1
   WRITE(Out_Unit,'(a,i2,a,ES16.7,"(Ha)",ES16.7,"(cm^-1)")')"Y(",i,")=",Yfit(i+1),Yfit(i+1)/cmm1toau
   WRITE(Spectra_Unit,'(a,i3,a,ES23.14,"(Ha)",ES23.14,"(cm^-1)")')"  i=",i,"  Yfit=",Yfit(i+1),Yfit(i+1)/cmm1toau
ENDDO
IF(jval==0)THEN
   !WRITE(Out_Unit,'(2A6,5A17,",",5A17)')"elec,","jval,","Yfit(0)(Ha),","Yfit(1)(Ha),","Yfit(2)(Ha),","Yfit(3)(Ha),","Yfit(4)(Ha),","Yfit(0)(cm^-1),","Yfit(1)(cm^-1),","Yfit(2)(cm^-1),","Yfit(3)(cm^-1),","Yfit(4)(cm^-1),"
   !WRITE(Spectra_Unit,'(2A6,5A17,",",5A17)')"elec,","jval,","Yfit(0)(Ha),","Yfit(1)(Ha),","Yfit(2)(Ha),","Yfit(3)(Ha),","Yfit(4)(Ha),","Yfit(0)(cm^-1),","Yfit(1)(cm^-1),","Yfit(2)(cm^-1),","Yfit(3)(cm^-1),","Yfit(4)(cm^-1),"
   WRITE(YfitPltDaf_Unit,'(2A6,5A17,",",5A17)')"elec,","jval,","Yfit(0)(Ha),","Yfit(1)(Ha),","Yfit(2)(Ha),","Yfit(3)(Ha),","Yfit(4)(Ha),","Yfit(0)(cm^-1),","Yfit(1)(cm^-1),","Yfit(2)(cm^-1),","Yfit(3)(cm^-1),","Yfit(4)(cm^-1),"
ENDIF
WRITE(YfitPltDaf_Unit,'(2(I5,","),5(ES16.7,","),",",5(ES16.7,","))')elec, jval,(Yfit(i+1),i=0,nYterms-1),(Yfit(i+1)/cmm1toau,i=0,nYterms-1)
WRITE(YfitBin_Unit)elec,jval,numnu,nbound,nquasi,nYTerms
WRITE(YfitBin_Unit)"DAF        ",DiatomicPot, Iarran, elec, jval, (Yfit(i+1),i=0,nYterms-1)

WRITE(DAFBin_Unit)(Yfit(i+1),i=0,nYterms-1)
FLUSH(DAFBin_Unit)
DEALLOCATE(Yfit,AMAT)

neig=nbound+nquasi
!CALL TimeInfoSub('BoundDAF','  End BoundDAF')

RETURN
END SUBROUTINE BoundDAF
