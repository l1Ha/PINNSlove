SUBROUTINE LogProp(jthenergy, working, dworking, ncall, smx_time, mcall, &
                    sge_time, ovrfile, refj, prty, nrho)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! modification of main_logder to treat photodissociation
!  M. Braunstein, Spring 1991
!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!  working is the file Rvwxyz in the working subdirectory into which the
!  R matrix will be put.
!  dworking is the file REvwxyz in the working subdirectory, etc.
!  version.
!
!  Changed number of Sfelevl_logder from elvl_unit to elev_unit.
!  Modified CALL to propph1() on 9/15/97 to conform to updates in that
!     routine: ipiv,ui,work now passed as arguments.
!<<trp
!
USE Narran_Module
USE FileUnits_MODULE
USE convrsns_MODULE
USE Masses_Module
USE prntng_MODULE
USE sublst_MODULE
USE csbasis_MODULE
USE TotalEng_Module
!      USE doenlvls_MODULE
USE dip_MODULE
USE engpro_MODULE
USE Opts_Module
USE InputFile_Module
USE pbasis_MODULE
USE Energy_Module
USE BasisRegions_Module
USE EDeriv_Module
USE LogDer_Module
USE APHchl_Module, ONLY: Neigmin_Sum, Neigmin_Save
!

IMPLICIT NONE
LOGICAL, Parameter :: Prt_Masses=.False.
REAL(Kind=WP_Kind) eredev, smx_time, sge_time, xktot, totmas
INTEGER nchanl, kkenergy, jtot, ncall, mcall, jthenergy, i, ii, iaph, IErr
CHARACTER(LEN=5) mlabel
CHARACTER(LEN= 8) today
CHARACTER(LEN=10) hour
CHARACTER(LEN=20) TypeOfMat, approx, proj, symm
CHARACTER(LEN=200) ovrfile
CHARACTER(LEN=*) working, dworking
LOGICAL jeven, nsymc(narran)
INTEGER lam_n(0:1000), lamfst, lamlst, ndim
INTEGER minmega, megacoup
INTEGER parity, megamax, jref, refj, jfact, prty, nrho
REAL(Kind=WP_Kind)  xtrafac, ipfac, tmursi, sqti, xlamp, xlampp, f_cor, f_asy
INTEGER, ALLOCATABLE :: chanl(:), elect(:), nvib(:), jrot(:)
INTEGER, ALLOCATABLE :: lorb(:), ipiv(:), tee(:), lambda(:)
REAL(Kind=WP_Kind), ALLOCATABLE :: R_Mat(:,:), w(:,:)
REAL(Kind=WP_Kind), ALLOCATABLE :: vecnew(:), tstore(:,:)
REAL(Kind=WP_Kind), ALLOCATABLE :: RE_Mat(:,:), zp(:,:)
REAL(Kind=WP_Kind), ALLOCATABLE :: w2(:,:),     zedz(:,:)
REAL(Kind=WP_Kind), ALLOCATABLE :: g0(:,:),     g1(:,:)
REAL(Kind=WP_Kind), ALLOCATABLE :: g2(:,:),     w0(:,:)
REAL(Kind=WP_Kind), ALLOCATABLE :: w1(:,:),     w2b(:,:)
REAL(Kind=WP_Kind), ALLOCATABLE :: diag(:,:),   cor(:,:)
REAL(Kind=WP_Kind), ALLOCATABLE :: asym(:,:),  dmu(:,:)
REAL(Kind=WP_Kind), ALLOCATABLE :: energies(:)
REAL(Kind=WP_Kind), ALLOCATABLE :: psi(:), gamma(:),  xksq(:)
REAL(Kind=WP_Kind), ALLOCATABLE :: zeta(:), omega(:), new(:)
REAL(Kind=WP_Kind), ALLOCATABLE :: ui(:,:), work(:)
!
COMMON /utotal/ xktot
!      COMMON / basis / nlam(0:1000), lammin, lammax, megamin, maxmega, coupmega
COMMON /totalj/ jfact
COMMON /xtra/ xtrafac, f_cor(0:1000), f_asy(0:1000)

NAMELIST/basis_n/ lam_n, lamfst, lamlst
NAMELIST/dipole2/ jf,prf,lamf,mf,jd,prd,lamd,md,ji,pri,lami, mi,nsfuncu,nsfuncl,lphoto,iaxis,enddip
NAMELIST/momentum/ jtot, parity, nsymc, jeven, jref, megacoup, megamax
EXTERNAL date, time, vsets, mxout, mxoutd

!----------------------------------------------------------------------- 
!   efirst   = Minimum Total Energy
!   deltaeng = Increment for the Total Energy
!   nenergy  = Number of Total Energies
!-----------------------------------------------------------------------
OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile,status='old')
   READ(In_Unit,NML=TotEnergy, IOSTAT=IERR)
   IF(IERR==-1)THEN
      CALL NameList_Default(Out_Unit, 'TotEnergy')
   ELSEIF(IERR/=0)THEN
      WRITE(Out_Unit,*)'ERROR with Namelist TotEnergy'
      WRITE(Msg_Unit,*)'ERROR with Namelist TotEnergy'
      STOP 'LogProp: TotEnergy'
   ENDIF
   IF(eV_Input)THEN
      Efirst=Efirst_eV/autoeV
      DeltaEng=DeltaEng_eV/autoeV
   ENDIF
   WRITE(Out_Unit,NML=TotEnergy)
CLOSE(Unit=In_Unit)
Approx="Not Set"
Proj="Not Set"
Symm="Not Set"

!-----------------------------------------------------------------------
! debug printing options
!
!     nsubs ...... number of routines which will have their debug
!                  WRITEs turned on
!
! The following variables must be defined for each routine
! for which you wish to turn on debug WRITEs:
!
!    subs(i) ...... name of the routine to turn on debug WRITEs
!                   must be a Hollerith CHARACTER string of no more than
!                   eight CHARACTERs
!    iprt(i) ...... what level of debug WRITEs:
!                    0=none
!                    1=minimum
!                    2=moderate
!                    3=maximum
!                    4=combination (must read in iter(j,i), j=1 to 3)
!   iter(j,i) ...... how many times to loop thru routine i with the
!                    jth debug level turned on (j=1 to 3)
!-----------------------------------------------------------------------
OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile,status='old')
   READ(In_Unit,NML=debug,IOSTAT=IERR)
   IF(IERR==-1)THEN
      option='selected'                   
      iprtall=  1
      nsubs=  100
      CALL NameList_Default(Out_Unit, 'Debug')
   ELSEIF(IERR/=0)THEN
      WRITE(Out_Unit,*)'ERROR with Namelist Debug'
      WRITE(Msg_Unit,*)'ERROR with Namelist Debug'
      STOP 'LogProp: Debug'
   ENDIF
   CALL DebugOptions(Out_Unit)   !WRITE(Out_Unit,NML=debug)
CLOSE(Unit=In_Unit)

!-----------------------------------------------------------------------
!  Main option parameters:
! LSFUNC    LOGICAL variable set equal to .true. if sfunc is to be
!           executed.
! LMATELM   LOGICAL variable set equal to .true. if matelm is to be
!           executed.
! LOVERLAP  LOGICAL variable set equal to .true. if overlap is to be
! CRAY      LOGICAL variable set equal to .true. if running on a Cray
! SCHEME    INTEGER variable: =1 => molecular problem (Gauss_Legendre qu
!                             =2 => coulomb problem (Gauss_Laguerre quad
! CSTEST    LOGICAL , usually false.  true tests basis for mega=1
! LMATCH    LOGICAL , is set equal to true to WRITE out information
!           for the FEM_to_ABM or DVR_to_ABM transformation.
!-----------------------------------------------------------------------
OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile,status='old')
   READ(In_Unit,  NML=options, IOSTAT=IERR)
   IF(IERR==-1)THEN
      CALL NameList_Default(Out_Unit, 'Options')
      scheme=1
      cstest=.False.
      lsfunc=.True.
      lmatelm=.True.
      loverlap=.True.
      lsample=.False.
      cray=.False.
      lmatch=.False.
      lplot=.False.
      lenlvls=.False. 
      lmesher=.True.
      laph3d=.False.
      sfuntype='ABM'
   ELSEIF(IERR/=0)THEN
      WRITE(Out_Unit,*)'ERROR with Namelist Options'
      WRITE(Msg_Unit,*)'ERROR with Namelist Options'
      STOP 'LogProp: Options'
   ENDIF
   WRITE(Out_Unit,NML=options)
CLOSE(Unit=In_Unit)

!-----------------------------------------------------------------------
! Read in the masses (carbon 12 mass units) for each atom.
!-----------------------------------------------------------------------
OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile,status='old')
   READ(In_Unit,  NML=AtomicMasses, IOSTAT=IERR)                            ! Atomic Masses
      IF(IERR==-1)THEN
         REWIND In_Unit
         READ(In_Unit,  NML=Atoms, IOSTAT=IERR)
         IF(Prt_Masses)THEN
         WRITE(Out_Unit,NML=Atoms)
         WRITE(Msg_Unit,NML=Atoms)
         WRITE(Out_Unit,*)
         ENDIF
         DO I=1,3
            CALL AtomicWeights(AtomicSymbol(I), MassNumber(I), AtomicNumber(I), Mass(I), Abundance(I), AtomicWeight(I), Notes(I))
            IF(Prt_Masses)THEN
            WRITE(Out_Unit, '(A,I3)')'Atomic Number = ', AtomicNumber(I)
            WRITE(Out_Unit, '(A,A3)')'Atomic Symbol = ', TRIM(AtomicSymbol(I))
            WRITE(Out_Unit, '(A,I3)')'Mass Number = '  , MassNumber(I)
            WRITE(Out_Unit, '(A,1PES23.15)')'Relative Atomic Mass = '  , Mass(I)
            WRITE(Out_Unit, '(A,1PES23.15)')'Isotopic Composition = '  , Abundance(I)
            WRITE(Out_Unit, '(A,1PES23.15)')'Standard Atomic Weight = ', AtomicWeight(I)
            WRITE(Out_Unit, '(A,A)')'Notes = ', TRIM(Notes(I))
            WRITE(Out_Unit,*)
            ENDIF
         ENDDO
      ELSEIF (IERR/=0)THEN
         WRITE(Msg_Unit,*)' Missing a necessary namelist: Atoms or AtomicMasses'
      ENDIF
   WRITE(Out_Unit,NML=AtomicMasses)
CLOSE(In_Unit)

!----------------------------------------------------------------------
! Specify the following:
!   JTOT   total angular momentum
!   IPAR   parity (IPAR=0=even-parity, IPAR=1=odd-parity)
!   NSYMC  is true only if the diatom in that channel is homonuclear.
!          NSYMC can be true for only one arrangement channel.
!  NPAR   true means use parity decoupling for the 3-body system,
!         false means no parity decoupling
!  JEVEN  true is even parity states of the asymptotic diatoms
!         false is odd parity states of the asymptotic diatoms
!----------------------------------------------------------------------
OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile,status='old')
   JRef=0
   MegaMax=0
   MegaCoup=0
   READ(In_Unit,  NML=momentum, IOSTAT=IERR)
   IF(IERR/=0)THEN
      WRITE(Out_Unit,*)'ERROR with Namelist momentum'
      WRITE(Msg_Unit,*)'ERROR with Namelist momentum'
      STOP 'LogProp: Momentum'
   ENDIF
   WRITE(Out_Unit,NML=momentum)
CLOSE(Unit=In_Unit,status='keep')

IF(megamax/=maxmega)THEN
    WRITE(Msg_Unit,*)' Warning: megamax/=maxmega.'
    WRITE(Msg_Unit,*)' megamax,maxmega=',megamax,maxmega
    STOP 'logprop'
ENDIF
IF(jref/=refj)THEN
    WRITE(Msg_Unit,*)' Warning: jref/=refj'
    WRITE(Msg_Unit,*)' jref,refj=',jref,refj
    STOP 'logprop'
ENDIF
IF(parity/=prty)THEN
    WRITE(Msg_Unit,*)' Warning: parity/=prty'
    WRITE(Msg_Unit,*)' parity,prty=',parity,prty
    STOP 'logprop'
ENDIF
IF(megacoup/=coupmega)THEN
    WRITE(Msg_Unit,*)' Warning: megacoup/=coupmega'
    WRITE(Msg_Unit,*)' megacoup,coupmega=',megacoup,coupmega
    STOP 'logprop'
ENDIF
IF(megacoup<megamax)THEN
    WRITE(Out_Unit,*)' Coupling ending at mega=', megacoup
ENDIF
!
minmega=0
IF((-1)**(jtot+parity)/=1) minmega=1

IF(minmega<megamin)THEN
    WRITE(Msg_Unit,*)' Warning: minmega<megamin'
    WRITE(Msg_Unit,*)' minmega,megamin=',minmega,megamin
    STOP 'logprop'
ENDIF

!----------------------------------------------------------------------
! LamMin = Minimum value for Lambda
! LamMax = Maximum value for Lambda
! Lam_N  = List Containing the Number of Channels for each Lambda
!----------------------------------------------------------------------
Lam_N=0
OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile,status='old')
   READ(In_Unit,NML=basis_n,IOSTAT=IERR)
   IF(IERR==-1)THEN
      WRITE(Out_Unit,*)'ERROR missing Namelist Basis_n'
      WRITE(Msg_Unit,*)'ERROR missing Namelist Basis_n'
      STOP 'LogProp: Basis_n'
   ELSEIF(IERR/=0)THEN
      WRITE(Out_Unit,*)'ERROR with Namelist Basis_n'
      WRITE(Msg_Unit,*)'ERROR with Namelist Basis_n'
      STOP 'LogProp: Basis_n'
   ENDIF
   WRITE(Out_Unit,NML=basis_n)
CLOSE(Unit=In_Unit)
IF (Neigmin_Sum == 0 .OR. Neigmin_Save(lamfst) == 0) THEN
   Neigmin_Sum = SUM(Lam_N)
   Neigmin_Save(lamfst:lamlst) = Lam_N(lamfst:lamlst)
END IF
!NChanl=SUM(Lam_N)
NChanl=Neigmin_Sum !tmpmodgregparker
! PORT (2026-10-02): in the two-electron-state (2x2) data path the
! propagation runs with 2N channels (block-diagonal 2N interaction matrix
! from potmatn_2x2).  The energy-independent storage still holds the
! single-state records, which potmatn_2x2 consumes with n = N.
BLOCK
  LOGICAL :: ps2x2_lp
  CALL Get2x2Switch(ps2x2_lp)
  IF(ps2x2_lp)THEN
     NChanl = 2*NChanl
     WRITE(Out_Unit,'(A,I6)') 'LogProp: 2x2 data path, NChanl doubled to ', NChanl
     WRITE(Msg_Unit,'(A,I6)') 'LogProp: 2x2 data path, NChanl doubled to ', NChanl
  ENDIF
END BLOCK
WRITE(Msg_Unit,'(A,I5,A,I5)')'NChanl=',Nchanl, ',  Total Angular Momentum (jtot)=',jtot
WRITE(Out_Unit,*)'NChanl=',Nchanl
lammin=lamfst
lammax=lamlst

ALLOCATE(chanl(NChanl), elect(NChanl), nvib(NChanl), jrot(NChanl))
ALLOCATE(lorb(NChanl), ipiv(NChanl), tee(NChanl), lambda(NChanl))
ALLOCATE(R_Mat(NChanl,NChanl), w(NChanl,NChanl))
ALLOCATE(vecnew(NChanl), tstore(NChanl,NChanl))
ALLOCATE(RE_Mat(NChanl,NChanl), zp(NChanl,NChanl))
ALLOCATE(w2(NChanl,NChanl),     zedz(NChanl,NChanl))
ALLOCATE(g0(NChanl,NChanl),     g1(NChanl,NChanl))
ALLOCATE(g2(NChanl,NChanl),     w0(NChanl,NChanl))
ALLOCATE(w1(NChanl,NChanl),     w2b(NChanl,NChanl))
ALLOCATE(diag(NChanl,NChanl),   cor(NChanl,NChanl))
ALLOCATE(asym(NChanl,NChanl),  dmu(NChanl,NChanl))
ALLOCATE(energies(NChanl))
ALLOCATE(psi(NChanl), gamma(NChanl),  xksq(NChanl))
ALLOCATE(zeta(NChanl), omega(NChanl), new(NChanl))
ALLOCATE(ui(NChanl,NChanl), work(NChanl))
!
IF(lammin<minmega) GOTO 240
IF(lammin>minmega)THEN
    WRITE(Msg_Unit,*)' Warning, lammin>minmega=',lammin,minmega
    WRITE(Out_Unit,*)' Warning, lammin>minmega=',lammin,minmega
ENDIF
IF(megacoup<lammax.AND.lammin/=lammax) GOTO 250
IF(lammax>megamax) GOTO 210
IF(lammax>jtot) GOTO 220
!
ndim=0
DO i=lammin,lammax
    nlam(i)=lam_n(i)
    ndim=max(nlam(i),ndim)
ENDDO
!Neigmin_sum=MAX(ndim,Neigmin_Sum)
IF(Neigmin_Sum<ndim)ndim=Neigmin_Sum !tmpmodgregparker
IF(ndim<=0) goto 260

chanl=0
elect=0
nvib=0
jrot=0
lorb=0
tee=0
lambda=0

iaph=0
DO i=lammin,lammax
    nlam(i)=Neigmin_Save(i)
    DO ii=1,nlam(i)
    iaph=iaph+1
    tee(iaph)=ii
    lambda(iaph)=i
    ENDDO
ENDDO

! PORT (2026-10-02): in the 2x2 data path nchanl = 2*iaph (the storage
! holds iaph single-state channels which potmatn_2x2 expands to 2N)
BLOCK
  LOGICAL :: ps2x2_lp2
  CALL Get2x2Switch(ps2x2_lp2)
  IF(ps2x2_lp2)THEN
     IF(iaph/=nchanl/2) goto 200
  ELSE
     IF(iaph/=nchanl) goto 200
  ENDIF
END BLOCK
IF(jthenergy==1)THEN
    OPEN(Unit=APH_Quant_Bin_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/APH_Quant.bin',form='unformatted', status='unknown')
        mlabel='aph  '
        ! PORT (2026-10-02): the stored quantum numbers keep the single-state
        ! count; the 2x2 data path only doubles the *propagation* channels
        BLOCK
          LOGICAL :: ps2x2_lp3
          INTEGER :: nqw
          CALL Get2x2Switch(ps2x2_lp3)
          nqw = nchanl
          IF(ps2x2_lp3) nqw = iaph
          WRITE(APH_Quant_Bin_Unit)mlabel
          WRITE(APH_Quant_Bin_Unit)jtot
          WRITE(APH_Quant_Bin_Unit)nqw
          WRITE(APH_Quant_Bin_Unit)(tee(iaph),iaph=1,nqw)
          WRITE(APH_Quant_Bin_Unit)(lambda(iaph),iaph=1,nqw)
        END BLOCK
    CLOSE(Unit=APH_Quant_Bin_Unit)
    IF(lenlvls)THEN
        OPEN(Unit=elev_unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'Sfelevl_logder')
    ENDIF
ENDIF


OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile,status='old')
   method='new_logder'
   READ(In_Unit,  NML=logder, IOSTAT=IERR)
   IF(IERR==-1)THEN
      CALL NameList_Default(Out_Unit, 'Logder')
      EDeriv=.True.
      nsteps=20
      iprtx='none  '
      init_rmat='Initialize  ' 
      rmat_file='R_matrix'
      method='new_logder'
   ELSEIF(IERR/=0)THEN
      WRITE(Out_Unit,*)'ERROR with Namelist Logder'
      WRITE(Msg_Unit,*)'ERROR with Namelist Logder'
      STOP 'LogProp: Logder'
   ENDIF
   WRITE(Out_Unit,NML=logder)
CLOSE(Unit=In_Unit)

IF(nrho<NSectors)nrho=NSectors
IF(nrho<NSectors) goto 230
!
! revised by X.L.
!      OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//'.dipole',status='old')
!      READ(In_Unit,NML=NML=dipole)
!      CLOSE(Unit=In_Unit)
!      WRITE(Out_Unit,NML=dipole)
lphoto=.false.
!

OPEN(Unit=work_unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//working, form='unformatted', status='unknown')

IF(EDeriv)THEN
   OPEN(Unit=dwork_unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//dworking, form='unformatted',status='unknown')
ENDIF

!
IF(lphoto)THEN
   OPEN(Unit=phmat_unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'R_Matrixph.bin', form='unformatted',status='unknown')
ENDIF
!
CALL constant

!-----------------------------------------------------------------------
! Compute the system reduced mass
!-----------------------------------------------------------------------
mass(1)=mass(1)*amutoau
mass(2)=mass(2)*amutoau
mass(3)=mass(3)*amutoau
totmas=mass(1)+mass(2)+mass(3)
usys=sqrt(mass(1)*mass(2)*mass(3)/totmas)
usys2=2.0d0*usys
!----------------------------------------------------------------------
! xtrafac: mega = 1 asymetric_top factor
!          = -1/(2mu)*(-1)**(jtot+parity)*jtot*(jtot+1)/4.d0
!
! f_cor:   coriolis factor
!          = -1/(2mu)[jtot(jtot+1) - mega(mega+1)]**(1/2)* mega =  0
!            *[1 + (-1)**(jtot+parity)]/2**(1/2)
!          = -1/(2mu)*[jtot*(jtot+1)-mega*(mega+1)]**(1/2) mega/=0
!
! f_asy:   asymetric_top factor
!          = 1/(2mu)*[jtot(jtot+1) - mega(mega+1)]**(1/2)* mega =  0
!            *[jtot(jtot+1) - (mega+1)(mega+2)]**(1/2)/4*
!            *[1 + (-1)**(jtot+parity)]/2**(1/2)
!          = 1/(2mu)*[jtot(jtot+1) - mega(mega+1)]**(1/2)*
!            *[jtot(jtot+1) - (mega+1)(mega+2)]**(1/2)/4   mega/=0
!----------------------------------------------------------------------
ipfac=(-1)**(jtot+parity)
tmursi= 1.0d0/usys2
sqti=1.0d0/sqrt(2.0d0)
!
xtrafac=-tmursi*jtot*(jtot+1)*ipfac/4.0d0
DO i=lammin,lammax
    f_cor(i)=0.0d0
    f_asy(i)=0.0d0
    IF(i<jtot)THEN
    xlamp=sqrt(jtot*(jtot+1.0d0)-i*(i+1.0d0))
    f_cor(i)=-xlamp*tmursi
    IF(i<(jtot-1))THEN
        xlampp=sqrt(jtot*(jtot+1.0d0)-(i+1.0d0)*(i+2.0d0))
        f_asy(i)=tmursi*xlamp*xlampp/4.0d0
    ENDIF
    ENDIF
ENDDO
f_cor(0)=f_cor(0)*(1+ipfac)*sqti
f_asy(0)=f_asy(0)*(1+ipfac)*sqti
!
jfact=jtot*(jtot+1)-jref*(jref+1)
!----------------------------------------------------------------------

!------------------------------------------------------------------
! Call system routines to determine the date and time.
!------------------------------------------------------------------

CALL dateh(today)
CALL timeh(hour)

!------------------------------------------------------------------
! Write the titles for this particular run.
!------------------------------------------------------------------

WRITE(Out_Unit,1131)
WRITE(Out_Unit,18)today
WRITE(Out_Unit,19)hour
!------------------------------------------------------------------
! The maximum energy will determine the wavelength
! of the most oscillatory state which will be used to set the
! step size in the log-derivative propagation.
!------------------------------------------------------------------
BLOCK
   CHARACTER(LEN=32) :: env_log
   INTEGER :: stat_log
   CALL get_environment_variable('APH3D_LOG_ENERGY', env_log, status=stat_log)
   IF (stat_log == 0 .AND. TRIM(ADJUSTL(env_log)) == '1') THEN
      engmax = eref + deltaeng
   ELSE
      engmax = efirst + (nenergy - 1) * deltaeng
   END IF
END BLOCK

!                   *****************************
!                   *   loop for each energy.   *
!                   *****************************
IF(lphoto)THEN
   OPEN(Unit=fph_unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'FirstEngph.bin',form='unformatted', status='old')
ENDIF
!      OPEN(Unit=frst_unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'Firsteng.bin',form='unformatted', status='old')
! revised by X.L.
!      CALL oli_read

OPEN(Unit=Ovr_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//ovrfile,form='unformatted', status='old')
REWIND(Ovr_Unit)
kkenergy=jthenergy
WRITE(Out_Unit,*)'kkenergy=',kkenergy, '     nenergy=',nenergy
!         READ(frst_unit) usys2
!------------------------------------------------------------------
! Determine the Total collision energy.
!------------------------------------------------------------------
BLOCK
   CHARACTER(LEN=32) :: env_log
   INTEGER :: stat_log
   REAL(Kind=WP_Kind) :: ecol_min, ecol_max, frac, e_coll
   CALL get_environment_variable('APH3D_LOG_ENERGY', env_log, status=stat_log)
   IF (stat_log == 0 .AND. TRIM(ADJUSTL(env_log)) == '1') THEN
      ecol_min = efirst
      ecol_max = deltaeng
      IF (nenergy > 1) THEN
         frac = DBLE(kkenergy - 1) / DBLE(nenergy - 1)
      ELSE
         frac = 0.0d0
      END IF
      e_coll = EXP(LOG(ecol_min) + frac * (LOG(ecol_max) - LOG(ecol_min)))
      Etot = eref + e_coll
   ELSE
      Etot = efirst + (kkenergy - 1) * deltaeng
   END IF
END BLOCK
xktot = usys2*Etot
eredev=Etot*autoev
WRITE(Out_Unit,'(A,ES23.15,A,ES23.15)')'energy (Ha) =',Etot, '   energy (eV) =',eredev
WRITE(Msg_Unit,'(A,ES21.13,A,ES21.13)')'energy (Ha) =',Etot, '   energy (eV) =',eredev
WRITE(Out_Unit,'(A,ES23.15,A,ES23.15)')'Etot=',Etot,'   xktot=',xktot
WRITE(Out_Unit,'(A,ES23.15)')' usys2=',usys2
!------------------------------------------------------------------
!  Set all elements of the LOG_derivative matrix initially to zero.
!  Then, set the diagonal elements to a large positive number.
!------------------------------------------------------------------
R_Mat=0.0d0
! PORT (2026-10-02): the hard-wall start value can be lowered with
! APH3D_WALL (default 1d30, the production value).  A smaller (still
! effectively infinite) wall greatly improves the conditioning of the
! first-sector propagation and is used to validate the 2N data path.
BLOCK
  REAL(Kind=WP_Kind) :: wallv
  CHARACTER(LEN=64) :: we
  INTEGER :: wl, ws
  CHARACTER(LEN=256) :: rmat_file_in
  INTEGER :: rmlen, rmstat, r_unit, nin, j
  REAL(Kind=WP_Kind) :: r_in
  LOGICAL :: rmat_exist
  wallv = 1d30
  CALL get_environment_variable('APH3D_WALL', we, wl, ws)
  IF(ws==0 .AND. wl>0) READ(we(1:wl),*) wallv
  IF(wallv /= 1d30) WRITE(Out_Unit,'(A,ES12.4)') 'LogProp: wall start value = ', wallv
  DO i=1,nchanl
     R_Mat(i,i)=wallv
  ENDDO
  ! PORT (2026-10-04, Route B DAT handoff): optionally read an initial R-matrix
  ! from an external file (e.g. from TwoState_Inner at r_DAT = 27.05 a0)
  rmat_file_in = ''
  CALL get_environment_variable('APH3D_RMAT_IN', rmat_file_in, rmlen, rmstat)
  IF(rmstat == 0 .AND. rmlen > 0) THEN
     INQUIRE(FILE=rmat_file_in(1:rmlen), EXIST=rmat_exist)
     IF(rmat_exist) THEN
        r_unit = 99
        OPEN(Unit=r_unit, File=rmat_file_in(1:rmlen), FORM='UNFORMATTED', STATUS='OLD')
        READ(r_unit) nin, r_in
        IF(nin == nchanl) THEN
           READ(r_unit) ((R_Mat(i,j), i=1,nchanl), j=1,nchanl)
           WRITE(Out_Unit,'(A,A,A,I5,A,F10.4)') 'LogProp: Successfully initialized R_Mat from ', &
                trim(rmat_file_in(1:rmlen)), ' (nin=', nin, ', r=', r_in, ')'
        ELSE
           WRITE(Out_Unit,'(A,I5,A,I5)') 'LogProp ERROR: APH3D_RMAT_IN dimension mismatch nin=', &
                nin, ' nchanl=', nchanl
           STOP 'LogProp: RMAT_IN dimension mismatch'
        ENDIF
        CLOSE(r_unit)
     ELSE
        WRITE(Out_Unit,*) 'LogProp WARNING: APH3D_RMAT_IN specified but file not found: ', &
             trim(rmat_file_in(1:rmlen))
     ENDIF
  ENDIF
END BLOCK
IF(EDeriv)THEN
   RE_Mat=0.0d0
   DO i=1,nchanl
      RE_Mat(i,i)=1d30
   ENDDO
ENDIF
!
! if we are photodissociating, initialize zeta matrix to zero
!
IF(lphoto)THEN
   zeta=0.d0
ENDIF
!--------------------------------------------------------------------
! Propagate the coupled equations in the APH region from start
! to finish if this region has nonzero length.
!--------------------------------------------------------------------
WRITE(Out_Unit,'("Using propph1 in logprop")')
CALL propph1(nchanl, R_Mat, w, vecnew, tstore, start,   &
             finish, nsteps, RE_Mat, zp, usys2, w2, zedz,   &
             eredev, EDeriv, iprtx, zeta,   &
             omega, new, g0, g1, g2, w0, w1, w2b,   &
             psi, dmu, energies, gamma,   &
             xksq, ncall, smx_time, mcall, sge_time,   &
             method, diag, cor, asym, ndim, ipiv, ui, work)
!--------------------------------------------------------------------
! Write the R_Matrix to a file subsequent for analysis.
!--------------------------------------------------------------------
WRITE(Out_Unit,*)'logprop writing R_Matrix'
TypeOfMat = 'R_Matrix'
CALL Write_Matrix2Old(R_Mat, energies, chanl, elect, nvib, jrot,   &
                  lorb, nchanl, Etot, kkenergy, TypeOfMat, jtot,   &
                  proj, approx, symm, work_unit)
CLOSE(Unit=work_unit)
!--------------------------------------------------------------------
! Write the derivative of R_Matrix with respect to E in a file
! for subsequent analysis.
!--------------------------------------------------------------------
IF(EDeriv)THEN
    TypeOfMat = 'RE_Matrix'
    CALL Write_Matrix2Old(RE_Mat, energies, chanl, elect, nvib, jrot,   &
                      lorb, nchanl, Etot, kkenergy, TypeOfMat, jtot,   &
                      proj, approx, symm, dwork_unit)
    CLOSE(Unit=dwork_unit)
ENDIF
!--------------------------------------------------------------------
! Write the gamma_2 matrix (in Manoloupolos' notation) to a file
!---------------------------------------------------------------------
IF(lphoto)THEN
WRITE(Out_Unit,*)'writing gamma_2 matrix to file'
WRITE(phmat_unit)(zeta(i),i=1,nchanl)
WRITE(Out_Unit,*)'zeta (1) = ', zeta (1), 'energy (eV) = ', eredev
ENDIF
!
iprtx='little'
IF(iprtx=='little'.or.iprtx=='medium')THEN
    WRITE(Out_Unit,*)'R_Matrix at ', finish,'Bohr', "  nchanl=", nchanl
    CALL MxOut(R_Mat, nchanl, nchanl)
    IF(EDeriv)THEN
          WRITE(Out_Unit,*)'Derivative of R_Matrix at', finish,'Bohr', "  nchanl=", nchanl
    CALL MxOut(RE_Mat, nchanl, nchanl)
    ENDIF
ELSEIF(iprtx=='full  ')THEN
    WRITE(Out_Unit,*)'R_Matrix at ', finish,'Bohr', "  nchanl=", nchanl
    CALL MxOutD(R_Mat, nchanl, nchanl, 0)
    IF(EDeriv)THEN
       WRITE(Out_Unit,*)'Derivative of R_Matrix at', finish,'Bohr', "  nchanl=", nchanl
       CALL MxOut(RE_Mat, nchanl, nchanl)
    ENDIF
ENDIF

!--------------------------------------------------------------------
! End of the energy do loop.
!--------------------------------------------------------------------
!      REWIND frst_unit
!
IF(lphoto)THEN
   REWIND fph_unit
ENDIF

CLOSE(Unit=Ovr_Unit,status='keep')
!
IF(lphoto)THEN
   CLOSE(Unit=fph_unit,status='keep')
ENDIF

DEALLOCATE(chanl, elect, nvib, jrot)
DEALLOCATE(lorb, ipiv, tee, lambda)
DEALLOCATE(R_Mat, w)
DEALLOCATE(vecnew, tstore)
DEALLOCATE(RE_Mat, zp)
DEALLOCATE(w2,   zedz)
DEALLOCATE(g0,     g1)
DEALLOCATE(g2,     w0)
DEALLOCATE(w1,    w2b)
DEALLOCATE(diag,  cor)
DEALLOCATE(asym,  dmu)
DEALLOCATE(energies)
DEALLOCATE(psi, gamma,  xksq)
DEALLOCATE(zeta, omega, new)
DEALLOCATE(ui, work)

!
19 FORMAT(1x,'Time: ',a10)
18 FORMAT(1x,'Date: ',a10)
1131 FORMAT(///)
1132 FORMAT(1x,a70)

RETURN
200  WRITE(Msg_Unit,*)' nchanl/=iaph of lam_n(i): nchanl, iaph=', nchanl, iaph
STOP 'logprop'
210  WRITE(Msg_Unit,*)' lammax>megamax: lammax,megamax=', lammax,megamax
STOP 'logprop'
220  WRITE(Msg_Unit,*)' lammax>jtot: lammax,jtot=', lammax,jtot
STOP 'logprop'
230  WRITE(Msg_Unit,*)' Error: nrho in abm < NSectors in .logder file'
WRITE(Msg_Unit,*)' nrho, NSectors=', nrho, NSectors
STOP 'logprop'
240  WRITE(Msg_Unit,*)' Error: lammin<minmega: lammin,minmega=', lammin, minmega
STOP 'logprop'
250  WRITE(Msg_Unit,*)' Error: megacoup<lammax.AND.lammin/=lammax'
WRITE(Msg_Unit,*)' megacoup, lammin, lammax= ', megacoup,lammin,lammax
STOP 'logprop'
260  WRITE(Out_Unit,*)' ndim le 0: ndim=', ndim
STOP 'logprop'
ENDSUBROUTINE LogProp
