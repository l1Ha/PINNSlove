SUBROUTINE readaph_ph1 (pfile, ofile)
!-----------------------------------------------------------------------
! Routine to read in all the input data
!-----------------------------------------------------------------------
! written by b.j. archer
!-----------------------------------------------------------------------
!
! updated Jan. 1991 by M. Braunstein to include photodissociation input
!

USE FileUnits_MODULE
USE nstate_MODULE
USE Narran_Module
USE parms_MODULE
USE dip_MODULE
USE Fbeta_MODULE
USE sfntyp_MODULE
USE crays_MODULE
USE elemnt_MODULE
USE Masses_Module
USE quantb_MODULE
USE GaussQuady_Module
USE pow_MODULE
USE Spectro_MODULE
USE Regins_Module
USE Title1_Module
USE aphchl_MODULE
USE opts_MODULE
USE totj_MODULE
USE PES_MODULE
USE prntng_MODULE
USE scalag_MODULE
USE Integrat_Module
USE qlam_MODULE
USE propag_MODULE
USE potzero_MODULE
USE csbasis_MODULE
USE engpro_MODULE
USE Storage_Module
USE sublst_MODULE
USE quad_MODULE
USE InputFile_MODULE
USE Energy_Module
USE Integrat_Module
USE Regins_Module
USE Boundary_Module
USE BasisRegions_Module
USE VIVS_Module
USE EDeriv_Module
USE Storage_Module
USE PotZero_Module
USE LogDer_Module
USE Stst_module
IMPLICIT NONE
LOGICAL, Parameter :: Prt_Masses=.False.
LOGICAL There
CHARACTER(LEN=*) pfile, ofile
CHARACTER(LEN=200) pmatfile, ovrfile
CHARACTER(LEN=200) sbfrst, jtotfile, corfile, asymfile
INTEGER i, ii, ithcall, iarran, IErr, vib
INTEGER lam_n(0:1000), lamfst, lamlst
EXTERNAL regons, popt, glegen, massfac
NAMELIST/quant/ minvib, maxvib, jmin, jmax
NAMELIST/basis_n/ lam_n, lamfst, lamlst
NAMELIST/ pmaovr/ pmatfile, ovrfile,sbfrst,jtotfile,corfile, asymfile
NAMELIST/ dipole2/ jf,prf,lamf,mf,jd,prd,lamd,md,ji,pri,lami, mi,nsfuncu,nsfuncl,lphoto,iaxis,enddip
!
!------------------------------------------------------------------
! Data statements
!------------------------------------------------------------------
!DATA little/.false./,medium/.false./,full/.false./, ithcall/0/,ithsub/0/
DATA ithcall /0/

npar = .true.
vzero = 0.d0
jeven = .true.
DO iarran = 1, narran
    chncof(iarran) = 1.d0
    betafrst(iarran) = 0.d0
ENDDO
cstest = .false.
lenlvls = .false.
sfuntype = 'ABM'
lphoto = .false.

!-----------------------------------------------------------------------
! Get the titles
!-----------------------------------------------------------------------
OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile,status='old')
   READ(In_Unit,  NML=Title_Labels, IOSTAT=IERR)
   IF(IERR==-1)THEN
      CALL NameList_Default(Out_Unit, 'Title_Labels')
      Title1='H+H2 reactive scattering using PK2 potential'
      Title2='channel(1)=HH+H,    channel(2)=HH+H     channel(3)=HH+H'
      Title3='H3 Porter-Karplus PES'
   ELSEIF(IERR/=0)THEN
      WRITE(Out_Unit,*)'ERROR with Namelist Title_Labels'
      WRITE(Msg_Unit,*)'ERROR with Namelist Title_Labels'
      STOP 'ReadAPH_Ph1: Title_Labels'
   ENDIF
CLOSE(Unit=In_Unit, status='keep')
WRITE(Out_Unit,'(A)')Title1
WRITE(Out_Unit,'(A)')Title2
WRITE(Out_Unit,'(A)')Title3
WRITE(Out_Unit,*)

!--------------------------------------------------------------------
! PES_Name is a label to uniquely identify the potential.
!--------------------------------------------------------------------
WRITE(Out_Unit,'(A)') PES_Name
!-----------------------------------------------------------------------
!  Main option parameters:
! CRAY      LOGICAL variable set equal to .true. if running on a Cray
! LENLVLS   LOGICAL variable set equal to .true. if the propagator
!           is to calculate the energy correlations
! SCHEME    INTEGER variable: =1 => molecular problem (Gauss_Legendre qu
!                             =2 => coulomb problem (Gauss_Laguerre quad
! SFUNTYPE  Hollerith string of length 3.  Method used to calculate
!           surface functions.  FEM=finite element method, DVR= discrete
!           variable representation, ABM=analytic basis method.
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
      STOP 'ReadAPH_Ph1: Options'
   ENDIF
   WRITE(Out_Unit,NML=options)
CLOSE(Unit=In_Unit, status='keep')

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
      STOP 'ReadAPH_Ph1: Debug'
   ENDIF
   CALL DebugOptions(Out_Unit)   !WRITE(Out_Unit,NML=debug)
CLOSE(Unit=In_Unit,status='keep')
!------------------------------------------------------------------
! Determine printing options.
!------------------------------------------------------------------
CALL popt('readin  ',little,medium,full,ithcall,ithsub)
little=.true.
IF(little)THEN
   WRITE(Out_Unit,options)
   WRITE(Out_Unit,*)'$DEBUG'
   WRITE(Out_Unit,*)"In routine ReadAPH_ph1"
   WRITE(Out_Unit,*)' NSUBS=',nsubs
   DO i=1,nsubs
      WRITE(Out_Unit,16)iprt(i), TRIM(subs(i))
16       FORMAT(" IPRT=",i4,3x,"SUBS=",A)
      IF(iprt(i)==4)THEN
         WRITE(Out_Unit,*)' ITER=',iter(1,i),iter(2,i),iter(3,i)
      ENDIF
   ENDDO
   WRITE(Out_Unit,*)'$end'
ENDIF
!-----------------------------------------------------------------------
! Read in the masses (carbon 12 mass units) for each atom.
! compute the system reduced mass
!-----------------------------------------------------------------------! Atomic Masses
INQUIRE (File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile, exist = there)
IF(there)THEN
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
ELSE
   WRITE(Out_Unit,  * ) 'Data file must exist: ', InputDIR(1:LEN(TRIM(InputDIR)))//InputFile
   WRITE(Msg_Unit,  * ) 'Data file must exist: ', InputDIR(1:LEN(TRIM(InputDIR)))//InputFile
   STOP "ReadAPH_Ph1"
ENDIF
CLOSE(Unit=In_Unit,status='keep')
CALL massfac(little)

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
   READ(In_Unit,  NML=momentum, IOSTAT=IERR)
   IF(IERR/=0)THEN
      WRITE(Out_Unit,*)'ERROR with Namelist momentum'
      WRITE(Msg_Unit,*)'ERROR with Namelist momentum'
      STOP 'ReadAPH_Ph1: Momentum'
   ENDIF
   WRITE(Out_Unit,NML=momentum)
CLOSE(Unit=In_Unit,status='keep')

!----------------------------------------------------------------------
!  Read in the propagation data:
! EFIRST   total energy for the first propagation
! DELTAENG energy increment
! NENERGY  number of total energies
!-----------------------------------------------------------------------
OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile,status='old')
   READ(In_Unit,  NML=TotEnergy, IOSTAT=IERR)
   IF(IERR==-1)THEN
      CALL NameList_Default(Out_Unit, 'TotEnergy')
   ELSEIF(IERR/=0)THEN
      WRITE(Out_Unit,*)'ERROR with Namelist TotEnergy'
      WRITE(Msg_Unit,*)'ERROR with Namelist TotEnergy'
      STOP 'ReadAllData: TotEnergy'
   ENDIF
   WRITE(Out_Unit,NML=TotEnergy)
CLOSE(Unit=In_Unit,status='keep')
emin = efirst
ejump = deltaeng
numengs = nenergy

!-----------------------------------------------------------------------
! If INTEGRAT(i) = true, then channel i is propagated.
!-----------------------------------------------------------------------
OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile,status='old')
   READ(In_Unit,  NML=IntChanl, IOSTAT=IERR)
   IF(IERR==-1)THEN
      CALL NameList_Default(Out_Unit, 'IntChanl')
      Integrat=[.True.,.True.,.True.]
   ELSEIF(IERR/=0)THEN
      WRITE(Out_Unit,*)'ERROR with Namelist IntChanl'
      WRITE(Msg_Unit,*)'ERROR with Namelist IntChanl'
      STOP 'ReadAPH_Ph1: IntChanl'
   ENDIF
   WRITE(Out_Unit,NML=IntChanl)
CLOSE(Unit=In_Unit,status='keep')
!-----------------------------------------------------------------------
! If IRDINDEP = true, read the energy independent propagation info from
!               disk on the first energy.  If false, do not read from
!               disk on the first energy.
! If IWRINDEP = true, WRITE the energy independent propagation info onto
!               disk on the first energy.  If false, do not WRITE
!               onto disk for the first energy.
!-----------------------------------------------------------------------
OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile,status='old')
   READ(In_Unit,  NML=storage, IOSTAT=IERR)
   IF(IERR==-1)THEN
      CALL NameList_Default(Out_Unit, 'Storage')
      irdindep=.false.
      iwrindep=.true.
   ELSEIF(IERR/=0)THEN
      WRITE(Out_Unit,*)'ERROR with Namelist Storage'
      WRITE(Msg_Unit,*)'ERROR with Namelist Storage'
      STOP 'ReadAPH_Ph1: Storage'
   ENDIF
   WRITE(Out_Unit,NML=storage)
CLOSE(Unit=In_Unit,status='keep')

!-----------------------------------------------------------------------
!  Read in data to determine the asymptotic vibrational-rotational
!  states of the isolated diatoms.
!  MINVIB   the minimum vibrational state for each arrangement channel.
!  MAXVIB   the maximum vibrational state for each arrangement channel.
!  JMIN     the minimum rotational state for each vibrational state
!           and each arrangement channel.
!  JMAX     the maximum rotational state for each vibrational state
!           and each arrangement channel.
!-----------------------------------------------------------------------
OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile,status='old')
   READ(In_Unit,  NML=quantum, IOSTAT=IERR)
   IF(IERR==-1)THEN
      CALL NameList_Default(Out_Unit, 'Quantum')
      minvib=[0,0,0]
      maxvib=[4,4,4]
      jmin=0
      jmax(0,1:3)=(Maxvib+1)*6 ! Use (maxvib+1)*(jmax)
   ELSEIF(IERR/=0)THEN
      WRITE(Out_Unit,*)'ERROR with Namelist Quantum'
      WRITE(Msg_Unit,*)'ERROR with Namelist Quantum'
      STOP 'ReadAPH_Ph1: Quantum'
   ENDIF
   DO iarran=1,narran
      nlegndre(iarran)=jmax(0,iarran)
      DO Vib=MinVib(iarran),MaxVib(iarran)
         jmax(Vib,iarran)=jmax(0,iarran)
      ENDDO
   ENDDO
   WRITE(Out_Unit,NML=quantum)
CLOSE(Unit=In_Unit,status='keep')
IF(little)THEN
   WRITE(Out_Unit,*)'$QUANTUM'
   DO i=1,narran
      IF(integrat(i))THEN
         WRITE(Out_Unit,*)' I=',i,' MINVIB=',minvib(i)
         WRITE(Out_Unit,*)' I=',i,' MAXVIB=',maxvib(i)
         WRITE(Out_Unit,*)' JMIN=',(jmin(ii,i,0), ii=minvib(i),maxvib(i))
         WRITE(Out_Unit,*)' JMAX=',(jmax(ii,i), ii=minvib(i),maxvib(i))
      ENDIF
   ENDDO
   WRITE(Out_Unit,*)'$end'
ENDIF

!-----------------------------------------------------------------------
!  Read in data necessary to generate the asymptotic basis and for
!  expanding the interaction potential so that interaction matrix
!  elements can be calculated analytically.
!  NOSCIL   number of harmonic oscillators basis functions.
!  NGLEGN   number of Gauss_Legendre points used to expand the
!           interaction potential in a Legendre polynomial.
!  NHERMT   numer of Gauss_Hermite points used in
!           calculating the Hamiltonian matrix elements.
!  RE       equilibrium position of the diatom for each arrangement
!           channels
!  RX       RX times RE is the position of the minimum  of the
!           parabola which defines the harmonic oscillator basis.
!           This parameter can usually be set equal to 1.1 for each
!           arrangement channel.  This parameter should be greater than
!           one because the long range part of the potential is softer
!           than short range part of the potential.
!  RALFA    RALFA times ALFA is curvature of the parabola which defines
!           the harmonic oscillator basis. This number should be
!           slightly less than one so that the harmonic oscillator basis
!           will have amplitude over a large enough range to adequately
!           expand the exact wavefunction.
!           This parameter can usually be set equal to 0.95 for each
!           arrangement channel.
!  NPOW     the order of the polynomial used to expand the vibrational
!           dependence of the interaction potential. The variable
!           x=(r-re)/re is used as the independent variable.
!-----------------------------------------------------------------------
OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile,status='old')
   READ(In_Unit,  NML=gauss, IOSTAT=IERR)
   IF(IERR==-1)THEN
      CALL NameList_Default(Out_Unit, 'Gauss')
      old_way= .True.
      npow=    [ 0, 0, 0]
      nglegn=  [46,46,46]
      nhermt=  [25,25,25]
      nlegndre=[ 6, 6, 6]
      intwt=   [ 1, 1, 1]
      rx=    [1.085000D0, 1.085000D0, 1.085000D0]
      re=    [1.401120D0, 1.401120D0, 1.401120D0]
      zeta=  [1.000000D0, 1.000000D0, 1.000000D0]
      delta= [0.010000D0, 0.010000D0, 0.010000D0]
      weau=  [2.00534D-02,2.00534D-02,2.00534D-02]
      wexeau=[5.52847D-04,5.52847D-04,5.52847D-04]
      ralpha=[3.166000D0, 3.166000D0, 3.166000D0]
      balpha=[1.000000D0, 1.000000D0, 1.000000D0]
      calpha=[0.920000D0, 0.920000D0, 0.920000D0]
      dalpha=[0.030000D0, 0.030000D0, 0.030000D0]
      anharm=[0.800000D0, 0.800000D0, 0.800000D0]
   ELSEIF(IERR/=0)THEN
      WRITE(Out_Unit,*)'ERROR with Namelist Gauss'
      WRITE(Msg_Unit,*)'ERROR with Namelist Gauss'
      STOP 'ReadAPH_Ph1: Gauss'
   ENDIF
   noscil=maxvib+1
   WRITE(Out_Unit,NML=gauss)
CLOSE(Unit=In_Unit,status='keep')

!-----------------------------------------------------------------------
!  Read in the spectroscopic constants to define the experimental
!  vibrational
!  rotational states of the isolated diatoms. These constants are
!  obtained from Herzberg's book Spectra of Diatomic Molecules and have
!  units of inverse centimeters. These are used to insure that one
!  has calculated accurate asymptotic wavefunctions.  The calculated
!  energies are used and not the 'experimental' energies.
!-----------------------------------------------------------------------
IF(scheme==1)THEN
!      READ(In_Unit,NML=spectro)
!      CLOSE(Unit=In_Unit,status='keep')
!      IF(little) WRITE(Out_Unit,spectro)
ENDIF
!-----------------------------------------------------------------------
! eta(i) modifies the argument of the Hylleras functions to vary how
!        long range the fit is in each channel, i.
! scalagr modifies the weights and zeros of the Gauss_Laguerre's
!         to allow good integration of all channels
!-----------------------------------------------------------------------
IF(scheme==2)THEN
   OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile,status='old')
   READ(In_Unit,  NML=etachanl, IOSTAT=IERR)
   IF(IERR==-1)THEN
      CALL NameList_Default(Out_Unit, 'EtaChanl')
   ELSEIF(IERR/=0)THEN
      WRITE(Out_Unit,*)'ERROR with Namelist EtaChanl'
      WRITE(Msg_Unit,*)'ERROR with Namelist EtaChanl'
      STOP 'ReadAPH_Ph1: EtaChanl'
   ENDIF
   WRITE(Out_Unit,NML=etachanl)
   CLOSE(Unit=In_Unit,status='keep')
ENDIF

!-----------------------------------------------------------------------
!  Read in the number of APH channels used in the propagation and the
!  values of lambda (projection of the total angular momentum on the
!  body-fixed z-axis) for each channel.
!-----------------------------------------------------------------------
OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile,status='old')
   Lam_N=0
   READ(In_Unit,NML=basis_n,IOSTAT=IERR)
   IF(IERR==-1)THEN
      WRITE(Out_Unit,*)'ERROR missing Namelist Basis_n'
      WRITE(Msg_Unit,*)'ERROR missing Namelist Basis_n'
      STOP 'ReadAPH_ph1: Basis_n'
   ELSEIF(IERR/=0)THEN
      WRITE(Out_Unit,*)'ERROR with Namelist Basis_n'
      WRITE(Msg_Unit,*)'ERROR with Namelist Basis_n'
      STOP 'ReadAPH_Ph1: Basis_n'
   ENDIF
   WRITE(Out_Unit,NML=basis_n)
   CLOSE(Unit=In_Unit)
   NAPH=SUM(Lam_N)
CLOSE(Unit=In_Unit,status='keep')
!-----------------------------------------------------------------------
!  Read in the propagation distances.
!  START    starting distance used in the propagation.
!  ENDAPH   distance which signals the end of the APH region.
!  ENDDELVE distance which signals the end of the Delves region.
!  SwitchToVIVS   distance to switch from the Log-derivative propagator to
!           the VIVS propagator.
!  FINISH   maximum propagation distance. Distance at which asymptotic
!          functions are used to obtain the S-matrix (scattering matrix)
!
!  NOTE the following inequalities must be satisfied:
!           START <= ENDAPH <= ENDDELVE <= FINISH
!  and
!           START <= SwitchToVIVS <= FINISH
!
!  The following schematic may be useful.
!
!
!              Log-derivative              VIVS
!          --------------------- -------------------------
!         |         |           |            |            |
!      start      endaph     SwitchToVIVS       enddelve      finish
!          --------- ------------------------ ------------
!             APH            Delves             Jacobi
!
!                       distance ----->
!-----------------------------------------------------------------------
OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile,status='old')
   READ(In_Unit,  NML=regins, IOSTAT=IERR)
   IF(IERR==-1)THEN
      CALL NameList_Default(Out_Unit, 'Regins')
      start=2.1d0
      SwitchToVIVS=8.1d0
      finish=8.1d0
      endaph=8.1d0
      enddelve=8.1d0
   ELSEIF(IERR/=0)THEN
      WRITE(Out_Unit,*)'ERROR with Namelist Regins'
      WRITE(Msg_Unit,*)'ERROR with Namelist Regins'
      STOP 'ReadAPH_Ph1: Regins'
   ENDIF
   WRITE(Out_Unit,NMl=Regins)
CLOSE(Unit=In_Unit,status='keep')

CALL regons
!-----------------------------------------------------------------------
!  Read in nsteps used to control the accuracy in the log-derivative
!  NSTEPS   number of steps per asymptotic wavelength at the highest
!           total energy. NSTEP=20 is usually adequate and NSTEP=40
!           usually gives very accurate S-matrix elements.
!           For energies less than the higest total energy the
!           scattering wavefunction will be calculated more accurately.
!-----------------------------------------------------------------------
    WRITE(Out_Unit,*)'start,SwitchToVIVS=',start,SwitchToVIVS
    IF(start<SwitchToVIVS)THEN
       OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile,status='old')
       READ(In_Unit,NML=pmaovr, IOSTAT=IERR)
   IF(IERR==-1)THEN
      CALL NameList_Default(Out_Unit, 'PmaOvr')
      pmatfile='BinOut/PotMatrx.bin'
      ovrfile='BinOut/Ovrlp.bin'
      sbfrst='BinOut/SetBasis_First.bin'
      jtotfile='BinOut/Diag_Jtot.bin'
      corfile='BinOut/Coriolis.bin'
      asymfile='BinOut/Asym_Top.bin'
   ELSEIF(IERR/=0)THEN
      WRITE(Out_Unit,*)'ERROR with Namelist PmaOvr'
      WRITE(Msg_Unit,*)'ERROR with Namelist PmaOvr'
      STOP 'ReadAllData: PmaOvr'
   ENDIF
       WRITE(Out_Unit,NML=pmaovr)
       CLOSE(Unit=In_Unit,status='keep')
       pfile = pmatfile
       ofile = ovrfile
       WRITE(Out_Unit,*)'Pmatrx info. will be read from: ',pfile
       WRITE(Out_Unit,*)'Ovrlp info. will be read from: ',ofile
       OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile,status='old')
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
      STOP 'ReadAPH_Ph1: Logder'
   ENDIF
          WRITE(Out_Unit,NML=logder)
       CLOSE(Unit=In_Unit,status='keep')
    ENDIF

!-----------------------------------------------------------------------
!  Read in parameters used in VIVS for control of the accuracy in the
!  propagation.
!  DRNOW   initial interval size.
!  TOFF    tolerance parameter which requires that the perturbation
!          corrections are less than this number.  This parameter
!          is roughly the accuracy of the resulting S-matrix elements
!          at the lowest energy.  S-matrix elements at higher energies
!          are probably calculated more accurately because the
!          perturbation corrections are smaller.
!  IAPLHA  number of step per interval. Usually set equal to 4.
!-----------------------------------------------------------------------
    IF(SwitchToVIVS<finish)THEN
       OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile,status='old')
       READ(In_Unit,  NML=regins, IOSTAT=IERR)
   IF(IERR==-1)THEN
      CALL NameList_Default(Out_Unit, 'Regins')
      start=2.1d0
      SwitchToVIVS=8.1d0
      finish=8.1d0
      endaph=8.1d0
      enddelve=8.1d0
   ELSEIF(IERR/=0)THEN
      WRITE(Out_Unit,*)'ERROR with Namelist Regins'
      WRITE(Msg_Unit,*)'ERROR with Namelist Regins'
      STOP 'ReadAPH_Ph1: Regins'
   ENDIF
          WRITE(Out_Unit,NML=vivs)
       CLOSE(Unit=In_Unit,status='keep')
    ENDIF

!-----------------------------------------------------------------------
! Read in parameters for photodissociation calculation.
! LPHOTO - if true then treat photodissociation
! JF - total final angular momentum
! MF - projection of final total final angular momentum (space frame)
! LAMF - projection of final total final angular momentum (body frame)
! PRF - parity of final D function
! JD - photon angular momentum (1)
! MD - projection of photon angular momentum (space frame)
! LAMD - projection of photon angular momentum (body frame)
! PRD - parity photon D function
! JI - total initial angular momentum
! MI - projection of initial total final angular momentum (space frame)
! LAMI - projection of initial total final angular momentum (body frame)
! PRI - parity of initial D function
! NSFUNCU - number of basis functions retained for upper state
! NSFUNCL - number of basis functions retained for lower state
!--------------------------------------------------------------------
lphoto=.false.
IF(lphoto)THEN
   OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile,status='old')
      READ(In_Unit,NML=dipole2,IOSTAT=IERR)
   IF(IERR==-1)THEN
      LPhoto=.False.
      CALL NameList_Default(Out_Unit, 'Dipole2')
   ELSEIF(IERR/=0)THEN
      WRITE(Out_Unit,*)'ERROR with Namelist Dipole2'
      WRITE(Msg_Unit,*)'ERROR with Namelist Dipole2'
      STOP 'ReadAPH_Ph1: Dipole2'
   ENDIF
      WRITE(Out_Unit,NML=dipole2)
   CLOSE(Unit=In_Unit)
ENDIF

RETURN
!--------------------------***end-readin***----------------------------
END SUBROUTINE ReadAPH_Ph1
