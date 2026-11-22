SUBROUTINE readaph
!
!-----------------------------------------------------------------------
! Routine to read in all the input data
!-----------------------------------------------------------------------
! writtenby b.j. archer
!
!this routine is called by:
!         aph_to_delves
!this routine calls
!         popt,regons,massfac
!-----------------------------------------------------------------------
USE fileunits_Module
USE Numeric_Kinds_Module
USE Narran_Module
USE Parms_Module
USE STST_Module
USE fbeta_Module
USE sfntyp_Module
USE crays_Module
USE elemnt_Module
USE doenlvls_Module
USE Masses_Module
USE quantb_Module
USE pow_Module
USE Spectro_Module
USE Regins_Module
USE Title1_Module
USE aphchl_Module
USE opts_Module
USE totj_Module
USE PES_MODULE
USE prntng_Module
USE scalag_Module
USE Integrat_Module
USE qlam_Module
USE propag_Module
USE potzero_Module
USE CSBasis_Module
USE engpro_Module
USE Storage_Module
USE sublst_Module
USE quad_Module
USE Energy_Module
USE VIVS_Module
USE LogDer_Module
USE GaussQuady_Module
USE FBR_Module
USE InputFile_Module
USE PES_MODULE
USE EDeriv_Module

IMPLICIT NONE
INTEGER i, ii, ithcall, iarran, ntt, nrot, vib, IERR
INTEGER lam_n(0:1000), lamfst, lamlst
DATA IthCall/0/

EXTERNAL regons, popt, glegen, massfac
NAMELIST/basis_n/ lam_n, lamfst, lamlst
!------------------------------------------------------------------
! Data statements
!------------------------------------------------------------------
!DATA little/.true./,medium/.false./,full/.false./, ithcall/0/,ithsub/0/
npar = .true.
vzero = 0.d0
jeven = .true.
DO  iarran = 1,narran
   chncof(iarran) = 1.d0
   betafrst(iarran) = 0.d0
ENDDO
cstest = .false.
!lenlvls = .false.
sfuntype = 'ABM'


!-----------------------------------------------------------------------
! Get the titles
!-----------------------------------------------------------------------
OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile,form='formatted',status='old')
READ(In_Unit,NML=Title_Labels,IOSTAT=IERR)
   IF(IERR==-1)THEN
      CALL NameList_Default(Out_Unit, 'Title_Labels')    
   ELSEIF(IERR/=0)THEN
      WRITE(Out_Unit,*)'ERROR with Namelist Title_Labels'
      WRITE(Msg_Unit,*)'ERROR with Namelist Title_Labels'
      STOP 'ReadABM: Title_Labels'
   ENDIF
CLOSE(Unit=In_Unit,status='keep')
WRITE(Out_Unit,'(A)')Title1
WRITE(Out_Unit,'(A)')Title2
WRITE(Out_Unit,'(A)')Title3
WRITE(Out_Unit,*)
!--------------------------------------------------------------------
! Read in a label to uniquely identify the potential. (25 CHARACTERs)
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
READ(In_Unit,NML=Options,IOSTAT=IERR)
IF(IERR==-1)THEN
      CALL NameList_Default(Out_Unit, 'Options')
   ELSEIF(IERR/=0)THEN
      WRITE(Out_Unit,*)'ERROR with Namelist Options'
      WRITE(Msg_Unit,*)'ERROR with Namelist Options'
      STOP 'ReadABM: Options'
   ENDIF
CLOSE(Unit=In_Unit,status='keep')
WRITE(Out_Unit,NML=Options)
!-----------------------------------------------------------------------
! debug printing Options
!
!     nsubs ...... number of routines which will have their debug
!                  WRITEs turned on
!
! The following variables must be defined for each routine
! for which you wish to turn on debug WRITEs:
!
!    subs(i) ...... name of the routine to turn on debug WRITEs;
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
OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile, form='formatted',status='old')
READ(In_Unit,NML=debug,IOSTAT=IERR)
IF(IERR==-1)THEN
      CALL NameList_Default(Out_Unit, 'Debug')
   ELSEIF(IERR/=0)THEN
      WRITE(Out_Unit,*)'ERROR with Namelist Debug'
      WRITE(Msg_Unit,*)'ERROR with Namelist Debug'
      STOP 'ReadABM: Debug'
   ENDIF
CLOSE(Unit=In_Unit,status='keep')
CALL DebugOptions(Out_Unit)   !WRITE(Out_Unit,NML=debug)

!------------------------------------------------------------------
! Determine printing Options.
!------------------------------------------------------------------
CALL popt('readin  ',little,medium,full,ithcall,ithsub)
IF(little)THEN
  WRITE(Out_Unit,Options)
  WRITE(Out_Unit,*)"In routine ReadAPH"
  WRITE(Out_Unit,*)'$DEBUG'
  WRITE(Out_Unit,*)' NSUBS=',nsubs
  DO i=1,nsubs
    WRITE(Out_Unit,'(" IPRT=",i4,3x,"SUBS=",A)')iprt(i), TRIM(subs(i))
    IF(iprt(i)==4)THEN
      WRITE(Out_Unit,*)' ITER=',iter(1,i),iter(2,i),iter(3,i)
    ENDIF
  ENDDO
  WRITE(Out_Unit,*)'$END'
ENDIF

!----------------------------------------------------------------------
! Specify the following:
!   JTOT   total angular momentum
!   PARITY   parity (PARITY=0=even-parity, PARITY=1=odd-parity)
!   NSYMC  is true only if the diatom in that channel is homonuclear.
!          NSYMC can be true for only one arrangement channel.
!  NPAR   true means use parity decoupling for the 3-body system,
!         false means no parity decoupling
!  JEVEN  true is even parity states of the asymptotic diatoms
!         false is odd parity states of the asymptotic diatoms
!----------------------------------------------------------------------
OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile, form='formatted',status='old')
READ(In_Unit,NML=momentum,IOSTAT=IERR)
IF(IERR==-1)THEN
      CALL NameList_Default(Out_Unit, 'Momentum')
   ELSEIF(IERR/=0)THEN
      WRITE(Out_Unit,*)'ERROR with Namelist Momentum'
      WRITE(Msg_Unit,*)'ERROR with Namelist Momentum'
      STOP 'ReadABM: Momentum'
   ENDIF
CLOSE(Unit=In_Unit,status='keep')
WRITE(Out_Unit,NML=momentum)
symmetry = nsymc(1)
mxmega=jtot
!----------------------------------------------------------------------
! Set the zero of the potential (in Hartree's), if needed.
!----------------------------------------------------------------------
OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile,form='formatted',status='old')
READ(In_Unit,NML=potzero,IOSTAT=IERR)
IF(IERR==-1)THEN
      CALL NameList_Default(Out_Unit, 'PotZero')
   ELSEIF(IERR/=0)THEN
      WRITE(Out_Unit,*)'ERROR with Namelist PotZero'
      WRITE(Msg_Unit,*)'ERROR with Namelist PotZero'
      STOP 'ReadABM: PotZero'
   ENDIF
WRITE(Out_Unit,NML=potzero)
CLOSE(Unit=In_Unit,status='keep')
WRITE(Out_Unit,NML=potzero)

!----------------------------------------------------------------------
!  Read in the propagation data:
! EFIRST   total energy for the first propagation
! DELTAENG energy increment
! NENERGY  number of total energies
!-----------------------------------------------------------------------
OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile,form='formatted',status='old')
CLOSE(Unit=In_Unit,status='keep')
emin = efirst
ejump = deltaeng
numengs = nenergy
IF(little) WRITE(Out_Unit,TotEnergy)
!-----------------------------------------------------------------------
! If INTEGRAT(i) = true, then channel i is propagated.
!-----------------------------------------------------------------------
OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile,form='formatted',status='old')
READ(In_Unit,NML=IntChanl,IOSTAT=IERR)
IF(IERR==-1)THEN
      CALL NameList_Default(Out_Unit, 'IntChanl')
   ELSEIF(IERR/=0)THEN
      WRITE(Out_Unit,*)'ERROR with Namelist IntChanl'
      WRITE(Msg_Unit,*)'ERROR with Namelist IntChanl'
      STOP 'ReadABM: IntChanl'
   ENDIF
CLOSE(Unit=In_Unit,status='keep')
WRITE(Out_Unit,NML=IntChanl)
!-----------------------------------------------------------------------
! If IRDINDEP = true, read the energy independent propagation info from
!               disk on the first energy.  If false, do not read from
!               disk on the first energy.
! If IWRINDEP = true, WRITE the energy independent propagation info onto
!               disk on the first energy.  If false, do not WRITE
!               onto disk for the first energy.
!-----------------------------------------------------------------------
OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile,form='formatted',status='old')
READ(In_Unit,NML=storage,IOSTAT=IERR)
IF(IERR==-1)THEN
      CALL NameList_Default(Out_Unit, 'Storage')
   ELSEIF(IERR/=0)THEN
      WRITE(Out_Unit,*)'ERROR with Namelist Storage'
      WRITE(Msg_Unit,*)'ERROR with Namelist Storage'
      STOP 'ReadABM: Storage'
   ENDIF
CLOSE(Unit=In_Unit,status='keep')
WRITE(Out_Unit,NML=storage)
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
OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile, form='formatted',status='old')
READ(In_Unit,NML=quantum,IOSTAT=IERR)
IF(IERR==-1)THEN
   CALL NameList_Default(Out_Unit, 'Quantum')
ELSEIF(IERR/=0)THEN
   WRITE(Out_Unit,*)'ERROR with Namelist Quantum'
   WRITE(Msg_Unit,*)'ERROR with Namelist Quantum'
   STOP 'ReadABM: Quantum'
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
  DO i=1,narran
  IF(integrat(i))THEN
    WRITE(Out_Unit,*)' I=',i,' MINVIB=',minvib(i)
    WRITE(Out_Unit,*)' I=',i,' MAXVIB=',maxvib(i)
    WRITE(Out_Unit,*)' JMIN=',(jmin(ii,i,0), ii=minvib(i),maxvib(i))
    WRITE(Out_Unit,*)' JMAX=',(jmax(ii,i), ii=minvib(i),maxvib(i))
  ENDIF
 ENDDO
  WRITE(Out_Unit,*)'$END'
ENDIF

MxBasis=0
DO iarran=1,narran
DO vib=MinVib(iarran),MaxVib(iarran)
   MxBasis=MxBasis+(jmax(vib,narran)+1)
ENDDO
ENDDO
MxBasiss=MxBasis
!----------------------------------------------------------------------
! IF symmetry decoupling is desired make sure quantum numbers are correct.
! First check vibrational quantum numbers.
!----------------------------------------------------------------------

IF(symmetry)THEN
   IF(minvib(2)/=minvib(3).or.maxvib(2)/=maxvib(3))THEN
      WRITE(Out_Unit,*)'The number of vibrational ', 'states must be ', 'equivalent for channels with identical atoms'
      STOP 'readaph'
   ENDIF

!----------------------------------------------------------------------
! Now check rotational quantum numbers.
!----------------------------------------------------------------------
   DO  vib = minvib(2), maxvib(2)
IF(jmin(vib,2,0)/=jmin(vib,3,0).or.jmax(vib,2)/=jmax(vib,3))THEN
   WRITE(Out_Unit,*)'The number of rotational states for ', &
       'each vibrational state must be equal ', &
       'for the channels with identical atoms'
   STOP 'readaph'
ENDIF
   ENDDO

   DO vib = minvib(1), maxvib(1)
IF(jeven)THEN
   IF(mod(jmin(vib,1,0),2)/=0 .or.mod(jmax(vib,1),2)/=0)THEN
      WRITE(Out_Unit,*)'jmin and jmax must be even for ',  'symmetry channel'
      STOP 'ReadABM'
   ENDIF
ELSE
   IF(mod(jmin(vib,1,0),2)/=1 .or.mod(jmax(vib,1),2)/=1)THEN
      WRITE(Out_Unit,*)'jmin and jmax must be odd for ', 'symmetry channel'
      STOP 'readaph'
   ENDIF
ENDIF
   ENDDO


MxBasis=0
MxBasiss=0
DO iarran=1,narran
   DO vib=MinVib(iarran),MaxVib(iarran)
      IF(iarran==1)THEN
         MxBasis =MxBasis +(jmax(vib,iarran)+1)
         MxBasiss=MxBasiss+(jmax(vib,iarran)/2+1)
      ELSEIF(iarran==3)THEN
         MxBasis=MxBasis+(jmax(vib,iarran)+1)
      ELSE
         MxBasis =MxBasis +(jmax(vib,iarran)+1)
         MxBasiss=MxBasiss+(jmax(vib,iarran)+1)
      ENDIF
   ENDDO
ENDDO
ENDIF

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

OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile, form='formatted',status='old')
READ(In_Unit,NML=gauss,IOSTAT=IERR)
IF(IERR==-1)THEN
   CALL NameList_Default(Out_Unit, 'Gauss')
ELSEIF(IERR/=0)THEN
   WRITE(Out_Unit,*)'ERROR with Namelist Gauss'
   WRITE(Msg_Unit,*)'ERROR with Namelist Gauss'
   STOP 'ReadABM: Gauss'
ENDIF 
noscil=maxvib+1
WRITE(Out_Unit,NML=gauss)
CLOSE(Unit=In_Unit,status='keep')
MxGlegn=MAXVAL(NGlegn)
MaxHermt=MAXVAL(NHermt)

!revised by X.L.
!maxhermt=max(maxhermt,21)
!MxGlegn=max(MxGlegn,51)
!  check size of hermite quadrature parms
  DO iarran=1,3
     IF(nhermt(iarran)>maxhermt)THEN
        WRITE(Out_Unit,*)'nhermt(iarran)>maxhermt'
        WRITE(Out_Unit,*)iarran,nhermt(iarran),maxhermt
        STOP 'readaph'
     ENDIF
 ENDDO

! need revised X.L.23-08-2006
    mxn=max(mxn,maxval(noscil))

!
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
READ(In_Unit,NML=Spectro,IOSTAT=IERR)
IF(IERR==-1)THEN
      CALL NameList_Default(Out_Unit, 'Spectro')
   ELSEIF(IERR/=0)THEN
      WRITE(Out_Unit,*)'ERROR with Namelist Spectro'
      WRITE(Msg_Unit,*)'ERROR with Namelist Spectro'
      STOP 'ReadABM: Spectro'
   ENDIF
CLOSE(Unit=In_Unit,status='keep')
WRITE(Out_Unit,NML=Spectro)
ENDIF
!-----------------------------------------------------------------------
! eta(i) modifies the argument of the Hylleras functions to vary how
!        long range the fit is in each channel, i.
! scalagr modifies the weights and zeros of the Gauss_Laguerre's
!         to allow good integration of all channels
!-----------------------------------------------------------------------
IF(scheme==2)THEN
   OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile,form='formatted',status='old')
   READ(In_Unit,NML=etachanl,IOSTAT=IERR)
   IF(IERR==-1)THEN
      CALL NameList_Default(Out_Unit, 'EtaChanl')
   ELSEIF(IERR/=0)THEN
      WRITE(Out_Unit,*)'ERROR with Namelist EtaChanl'
      WRITE(Msg_Unit,*)'ERROR with Namelist EtaChanl'
      STOP 'ReadABM: EtaChanl'
   ENDIF
   CLOSE(Unit=In_Unit,status='keep')
   WRITE(Out_Unit,NML=etachanl)
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
      STOP 'ReadAPH: Basis_n'
   ELSEIF(IERR/=0)THEN
      WRITE(Out_Unit,*)'ERROR with Namelist Basis_n'
      WRITE(Msg_Unit,*)'ERROR with Namelist Basis_n'
      STOP 'ReadAPH: Basis_n'
   ENDIF
   WRITE(Out_Unit,NML=basis_n)
   CLOSE(Unit=In_Unit)
   NAPH=SUM(Lam_N)
CLOSE(Unit=In_Unit,status='keep')
!----------------------------------------------------------------------
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
READ(In_Unit,NML=regins,IOSTAT=IERR)
IF(IERR==-1)THEN
      CALL NameList_Default(Out_Unit, 'Regins')
   ELSEIF(IERR/=0)THEN
      WRITE(Out_Unit,*)'ERROR with Namelist Regins'
      WRITE(Msg_Unit,*)'ERROR with Namelist Regins'
      STOP 'ReadABM: Regins'
   ENDIF
CLOSE(Unit=In_Unit,status='keep')
WRITE(Out_Unit,NML=regins)
CALL regons
!-----------------------------------------------------------------------
!  Read in nsteps used to control the accuracy in the log-derivative
!  NSTEPS   number of steps per asymptotic wavelength at the highest
!           total energy. NSTEP=20 is usually adequate and NSTEP=40
!           usually gives very accurate S-matrix elements.
!           For energies less than the higest total energy the
!           scattering wavefunction will be calculated more accurately.
!-----------------------------------------------------------------------
!BEGIN TEMPMOD
!IF(start<SwitchToVIVS)THEN
!   OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile,form='formatted',status='old')
!   READ(In_Unit,NML=logder)
!   CLOSE(Unit=In_Unit,status='keep')
!   WRITE(Out_Unit,NML=logder)
!ENDIF
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
!IF(SwitchToVIVS<finish)THEN
!   OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile,form='formatted',status='old')
!   READ(In_Unit,NML=vivs)
!   CLOSE(Unit=In_Unit,status='keep')
!   WRITE(Out_Unit,NML=vivs)
!ENDIF
!END TEMPMOD

!-----------------------------------------------------------------------
!  Read in parameters used to construct the fbr wavefunction in the
!  aph to delves transformation (IF DVR was used).
!  LMAX    maximum order of Legendre polynomial in theta expansion
!  MMAX    maximum order of trignometric function in chi expansion
!  NSFUNC  the number of surface functions written to the
!         sfunfbr file by the dvr code
!  NSYM    the symmetry of the trignometric function in chi
!          used in the dvr calculation
!-----------------------------------------------------------------------
   IF(sfuntype=='DVR')THEN
      OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile,form='formatted',status='old')
      READ(In_Unit,NML=fbr,IOSTAT=IERR)
   IF(IERR==-1)THEN
      NSfunc=105    !TNPMODGregParker
      lmax=20
      mmax=40
      nsym=1
      CALL NameList_Default(Out_Unit, 'FBR')
   ELSEIF(IERR/=0)THEN
      WRITE(Out_Unit,*)'ERROR with Namelist FBR'
      WRITE(Msg_Unit,*)'ERROR with Namelist FBR'
      STOP 'ReadAPH: FBR'
   ENDIF
      CLOSE(Unit=In_Unit,status='keep')
      WRITE(Out_Unit,NML=fbr)
!         nthetmax=lmax+1                          !Question: Is this already set?
!         nchimax=mmax                             !Question: Is this already set?
!         ndvptmax=nthetmax*nchimax                !Question: Is this already set?
    ENDIF
ntt=0
DO i=1,narran
   nrot=jmax(0,i) + 1
   IF(nsymc(i))nrot = (nrot+1)/2
   ntt=ntt + nrot*noscil(i)
ENDDO
nvibrot=ntt
nvbrthrt=nvibrot*maxhermt
maxosc=maxhermt
IF(ntt>nvibrot)THEN
   WRITE(Out_Unit,*) 'ntt=', ntt, ' is greater than nvibrot=', nvibrot
   WRITE(Out_Unit,*) ' increase nvibrot.parms'
   WRITE(Msg_Unit,*) 'ntt=', ntt, ' is greater than nvibrot=', nvibrot
   WRITE(Msg_Unit,*) ' increase nvibrot.parms'
   STOP 'readaph'
ENDIF
IF(little) WRITE(Out_Unit,*) 'ntt=', ntt, ' nvibrot=', nvibrot
mxang=0
DO i=1,narran
IF(i<3)THEN
IF(nglegn(i)<nglegn(i+1))THEN
mxang=mxang+nhermt(i)*nglegn(i+1)
ELSE
mxang=mxang+nhermt(i)*nglegn(i)
ENDIF
ELSE
   mxang=mxang+nhermt(i)*nglegn(i)
ENDIF
ENDDO
mxang=max(mxang,2016)
RETURN
!--------------------------***END-readin***----------------------------
ENDSUBROUTINE ReadAPH
