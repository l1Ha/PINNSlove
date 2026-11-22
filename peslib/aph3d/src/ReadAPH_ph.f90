SUBROUTINE readaph_ph
  USE Numeric_Kinds_Module
!
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
  USE STST_Module
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
  USE Energy_MODULE
  USE Regins_Module
  USE PotZero_Module
  IMPLICIT NONE
  LOGICAL, Parameter :: Prt_Masses=.False.
  LOGICAL There
  INTEGER  i, ii, ithcall, IErr
  INTEGER ntt,nrot,vib,iarran
!     >  alfecm, becm, decm, wecm, wexecm, weyecm,
!     >, reau(narran), vshift(3)

  !LOGICAL little, medium, full

  EXTERNAL regons, popt, glegen, massfac
  NAMELIST/dipole2/ jf,prf,lamf,mf,jd,prd,lamd,md,ji,pri,lami, &
                       mi,nsfuncu,nsfuncl,lphoto,iaxis,enddip
!
!------------------------------------------------------------------
! Data statements
!------------------------------------------------------------------
!      DATA little/.false./,medium/.false./,full/.false./, ithcall/0/,ithsub/0/

  npar = .true.
  vzero = 0.d0
  jeven = .true.
  chncof(1)=1.d0
  chncof(2)=1.d0
  chncof(3)=1.d0
  betafrst(1)=0.d0
  betafrst(2)=0.d0
  betafrst(3)=0.d0
  cstest = .false.
  lenlvls = .false.
  sfuntype = 'ABM'
  lphoto = .false.
!-----------------------------------------------------------------------
! Get the titles
!-----------------------------------------------------------------------
OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile,form='formatted',status='old')
   READ(In_Unit,NML=Title_Labels)
CLOSE(In_Unit)

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
  OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile,form='formatted',status='old')
  READ(In_Unit,NML=options, IOSTAT=IERR)
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
      STOP 'ReadAPH_Ph: Options'
   ENDIF
  CLOSE(Unit=In_Unit, status='keep')
  WRITE(Out_Unit,NML=options)
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
   OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile,form='formatted',status='old')
!      OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//'.debug_options',
!     >      form='formatted',status='old')
  READ(In_Unit,NML=debug, IOSTAT=IERR)
   IF(IERR==-1)THEN
      option='selected'                   
      iprtall=  1
      nsubs=  100
      CALL NameList_Default(Out_Unit, 'Debug')
   ELSEIF(IERR/=0)THEN
      WRITE(Out_Unit,*)'ERROR with Namelist Debug'
      WRITE(Msg_Unit,*)'ERROR with Namelist Debug'
      STOP 'ReadAPH_ph: Debug'
   ENDIF
  CLOSE(Unit=In_Unit,status='keep')
  CALL DebugOptions(Out_Unit)   !WRITE(Out_Unit,NML=debug)
!------------------------------------------------------------------
! Determine printing options.
!------------------------------------------------------------------
  CALL popt('readin  ',little,medium,full,ithcall,ithsub)
  IF(little)THEN
    WRITE(Out_Unit,options)
    WRITE(Out_Unit,*)'$DEBUG'
  WRITE(Out_Unit,*)"In routine ReadAPH_ph"
    WRITE(Out_Unit,*)' NSUBS=',nsubs
    DO 10 i=1,nsubs
      WRITE(Out_Unit,16)iprt(i), TRIM(subs(i))
16       FORMAT(" IPRT=",i4,3x,"SUBS=",A)
      IF(iprt(i)==4)THEN
        WRITE(Out_Unit,*)' ITER=',iter(1,i),iter(2,i),iter(3,i)
      ENDIF
10     CONTINUE
    WRITE(Out_Unit,*)'$end'
  ENDIF
!-----------------------------------------------------------------------
! Read in the masses (carbon 12 mass units) for each atom.
! compute the system reduced mass
!-----------------------------------------------------------------------
 INQUIRE (File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile, exist = there)
  IF(there)THEN
   OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile)
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
      CLOSE(In_Unit)
ELSE
   WRITE(Out_Unit,  * ) 'Data file must exist: ', InputDIR(1:LEN(TRIM(InputDIR)))//InputFile
   WRITE(Msg_Unit,  * ) 'Data file must exist: ', InputDIR(1:LEN(TRIM(InputDIR)))//InputFile
   STOP "ReadAllData"
ENDIF
  WRITE(Out_Unit,NML=AtomicMasses)
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
   OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile,form='formatted',status='old')
!      OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//'.define_system',
!     >      form='formatted',status='old')
  READ(In_Unit,NML=momentum, IOSTAT=IERR)
   IF(IERR/=0)THEN
      WRITE(Out_Unit,*)'ERROR with Namelist momentum'
      WRITE(Msg_Unit,*)'ERROR with Namelist momentum'
      STOP 'ReadAPH_Ph: Momentum'
   ENDIF
  CLOSE(Unit=In_Unit,status='keep')
  WRITE(Out_Unit,NML=momentum)
!----------------------------------------------------------------------
!  Read in the propagation data:
! EFIRST   total energy for the first propagation
! DELTAENG energy increment
! NENERGY  number of total energies
!-----------------------------------------------------------------------
  OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile,form='formatted',status='old')
  READ(In_Unit,NML=TotEnergy, IOSTAT=IERR)
   IF(IERR==-1)THEN
      CALL NameList_Default(Out_Unit, 'TotEnergy')
   ELSEIF(IERR/=0)THEN
      WRITE(Out_Unit,*)'ERROR with Namelist TotEnergy'
      WRITE(Msg_Unit,*)'ERROR with Namelist TotEnergy'
      STOP 'ReadAPH_Ph: TotEnergy'
   ENDIF
  CLOSE(Unit=In_Unit,status='keep')
  WRITE(Out_Unit,NML=TotEnergy)
  emin = efirst
  ejump = deltaeng
  numengs = nenergy
!-----------------------------------------------------------------------
! If INTEGRAT(i) = true, then channel i is propagated.
!-----------------------------------------------------------------------
  OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile,form='formatted',status='old')
  READ(In_Unit,NML=IntChanl, IOSTAT=IERR)
   IF(IERR==-1)THEN
      CALL NameList_Default(Out_Unit, 'IntChanl')
      Integrat=[.True.,.True.,.True.]
   ELSEIF(IERR/=0)THEN
      WRITE(Out_Unit,*)'ERROR with Namelist IntChanl'
      WRITE(Msg_Unit,*)'ERROR with Namelist IntChanl'
      STOP 'ReadAPH_Ph: IntChanl'
   ENDIF
  CLOSE(Unit=In_Unit,status='keep')
  IF(little) WRITE(Out_Unit,IntChanl)
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
   OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile,form='formatted',status='old')
!      OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//'.jacobi_basis',
!     >      form='formatted',status='old')
  READ(In_Unit,NML=quantum, IOSTAT=IERR)
   IF(IERR==-1)THEN
      CALL NameList_Default(Out_Unit, 'Quantum')
      minvib=[0,0,0]
      maxvib=[4,4,4]
      jmin=0
      jmax(0,1:3)=(Maxvib+[1,1,1])*6 ! Use (maxvib+1)*(jmax)
   ELSEIF(IERR/=0)THEN
      WRITE(Out_Unit,*)'ERROR with Namelist Quantum'
      WRITE(Msg_Unit,*)'ERROR with Namelist Quantum'
      STOP 'ReadAPH_Ph: Quantum'
   ENDIF
  CLOSE(Unit=In_Unit,status='keep')
  DO iarran=1,narran
     nlegndre(iarran)=jmax(0,iarran)
     DO Vib=MinVib(iarran),MaxVib(iarran)
        jmax(Vib,iarran)=jmax(0,iarran)
     ENDDO
  ENDDO
  WRITE(Out_Unit,NML=quantum)
  IF(little)THEN
    WRITE(Out_Unit,*)'$QUANTUM'
    DO 24 i=1,narran
    IF(integrat(i))THEN
!          WRITE(Msg_Unit,*)i,minvib(i)
!          WRITE(Msg_Unit,*)i,maxvib(i)
      WRITE(Out_Unit,*)' I=',i,' MINVIB=',minvib(i)
      WRITE(Out_Unit,*)' I=',i,' MAXVIB=',maxvib(i)
      WRITE(Out_Unit,*)' JMIN=',(jmin(ii,i,0), ii=minvib(i),maxvib(i))
      WRITE(Out_Unit,*)' JMAX=',(jmax(ii,i), ii=minvib(i),maxvib(i))
    ENDIF
24     CONTINUE
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
OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile,form='formatted',status='old')
!      OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//'.jacobi_basis',
!     >      form='formatted',status='old')
   READ(In_Unit,NML=gauss, IOSTAT=IERR)
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
      STOP 'ReadAPH_Ph: Gauss'
   ENDIF
   noscil=maxvib+1
   WRITE(Out_Unit,NML=gauss)
CLOSE(Unit=In_Unit,status='keep')
  ntt=0
  DO i=1,narran
     nrot=jmax(0,i) + 1
     IF(nsymc(i))nrot = (nrot+1)/2
     ntt=ntt + nrot*noscil(i)
  ENDDO
  nvibrot=ntt
  MaxHermt=MaxVal(nhermt)
  nvbrthrt=nvibrot*maxhermt

  mxl=0
  mxn=0
  DO iarran=1,narran
     DO vib=MinVib(iarran),MaxVib(iarran)
      IF(Mxl<jmax(vib,iarran))mxl=jmax(vib,iarran)
     ENDDO
     IF(mxn<maxvib(iarran))mxn=maxvib(iarran)
  ENDDO
  mxl=mxl+1
  mxn=mxn+1
! revisd by X.L. to make noscil>maxvib
  mxn=max(mxn,maxval(noscil))
!      STOP

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
   OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile,form='formatted',status='old')
!      OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//'.spectro',form='formatted',status='old')
  READ(In_Unit,NML=spectro, IOSTAT=IERR)
   IF(IERR==-1)THEN
      CALL NameList_Default(Out_Unit, 'Spectro')
      diatom=['H-H    ','H-H    ','H-H    ']
      wecm=  [0.4401213d+04, 0.4401213d+04, 0.4401213d+04]
      wexecm=[0.1213360d+03, 0.1213360d+03, 0.1213360d+03]
      weyecm=[0.8129000d+00, 0.8129000d+00, 0.8129000d+00]
      becm=  [0.6085300d+02, 0.6085300d+02, 0.6085300d+02]
      decm=  [0.4710000d-01, 0.4710000d-01, 0.4710000d-01]
      alfecm=[0.3062200d+01, 0.3062200d+00, 0.3062200d+00]
      reau=  [0.1401081d+01, 0.1401081d+01, 0.1401081d+01]
      vshift=[0.0000000d+00, 0.0000000d+00, 0.0000000d+00]
   ELSEIF(IERR/=0)THEN
      WRITE(Out_Unit,*)'ERROR with Namelist Spectro'
      WRITE(Msg_Unit,*)'ERROR with Namelist Spectro'
      STOP 'ReadAPH_Ph: Spectro'
   ENDIF
  CLOSE(Unit=In_Unit,status='keep')
  WRITE(Out_Unit,NML=spectro)
  ENDIF
!-----------------------------------------------------------------------
! eta(i) modifies the argument of the Hylleras functions to vary how
!        long range the fit is in each channel, i.
! scalagr modifies the weights and zeros of the Gauss_Laguerre's
!         to allow good integration of all channels
!-----------------------------------------------------------------------
     IF(scheme==2)THEN
         OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//'.etachanl',form='formatted',status='old')
         READ(In_Unit,NML=etachanl, IOSTAT=IERR)
   IF(IERR==-1)THEN
      CALL NameList_Default(Out_Unit, 'EtaChanl')
   ELSEIF(IERR/=0)THEN
      WRITE(Out_Unit,*)'ERROR with Namelist EtaChanl'
      WRITE(Msg_Unit,*)'ERROR with Namelist EtaChanl'
      STOP 'ReadAPH_Ph: EtaChanl'
   ENDIF
         CLOSE(Unit=In_Unit,status='keep')
         WRITE(Out_Unit,NML=etachanl)
      ENDIF
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
  OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile,form='formatted',status='old')
    READ(In_Unit,NML=Regins, IOSTAT=IERR)
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
      STOP 'ReadAPH_Ph: Regins'
   ENDIF
    WRITE(Out_Unit,NML=Regins)
  CLOSE(Unit=In_Unit,status='keep')

  CALL regons
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
     OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile,form='formatted',status='old')
!      OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//'.dipole',form='formatted',status='old')
     READ(In_Unit,NML=dipole2,IOSTAT=IERR)
   IF(IERR==-1)THEN
      LPhoto=.False.
      CALL NameList_Default(Out_Unit, 'Dipole2')
   ELSEIF(IERR/=0)THEN
      WRITE(Out_Unit,*)'ERROR with Namelist Dipole2'
      WRITE(Msg_Unit,*)'ERROR with Namelist Dipole2'
      STOP 'ReadAPH_Ph: Dipole2'
   ENDIF
     CLOSE(Unit=In_Unit,status='keep')
     WRITE(Out_Unit,NML=dipole2)
  ENDIF
  
  RETURN
!--------------------------***end-readin***----------------------------
  END SUBROUTINE ReadAPH_Ph
