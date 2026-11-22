SUBROUTINE ReadABM(numrho)
!
!            P U R P O S E   O F    S U B R O U T I N E
!   This routine reads in DATA from UNIT=In_Unit

!            O U T P U T    A R G U M E N T S
!  rhomin   minimum hyperradius
!  rhomax   maximum hyperradius
!  jtot     total angular momentum
!  parity   parity
!  megamin  minimum projection of the total angular momentum.
!  megamax  maximum projection of the total angular momentum.
!  naph     naph number of surface function desired.
!  lsfunc   true implies calculate surface functions.
!  loverlap true implies calculate overlap matrix elements.
!  lmatelm  true implies calculate coupling matrix elements.

! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

USE Parms_Module
USE TotJ_Module
USE FileUnits_Module
USE PES_MODULE
USE Quantb_Module
USE Masses_Module
USE QDiv_Module
USE QCase_Module
USE Converge_Module
USE Numbers_Module
USE SfnDis_Module
USE Crays_Module
USE Opts_Module
USE Gaussb_Module
USE Prntng_Module
USE SubLst_Module
USE Boundary_Module
USE InputFile_Module
USE BasisRegions_Module
USE Title1_Module
IMPLICIT NONE
LOGICAL There
INTEGER i, vib, mod, numrho, idiv, iarran, nth, ithdiv, ichdiv, nang, nch, IERR
INTEGER NAngCh(NArran+1),maxtemp
INTEGER, ALLOCATABLE:: MxBasisMega(:), MxBasissMega(:)
EXTERNAL popt, massfac
INTRINSIC mod
!-----------------------------------------------------------------------
!  Obtain printing options.
!-----------------------------------------------------------------------
!LOGICAL little, medium, full
INTEGER ithcall
DATA ithcall /0/    !, ithsub/0/
!DATA little /.True./, medium/.false./, full/.false./

!-----------------------------------------------------------------------
! Get the titles
!-----------------------------------------------------------------------
old_way = .true.
OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile, form='formatted',status='old')
   READ(In_Unit,NML=Title_Labels)
CLOSE(Unit=In_Unit)
WRITE(Out_Unit,'(A)')Title1
WRITE(Out_Unit,'(A)')Title2
WRITE(Out_Unit,'(A)')Title3
WRITE(Out_Unit,*)


!-----------------------------------------------------------------------
!  Main option parameters:
! LSFUNC    LOGICAL variable set equal to .true. IF sfunc is to be
!           executed.
! LMATELM   LOGICAL variable set equal to .true. IF matelm is to be
!           executed.
! LOVERLAP  LOGICAL variable set equal to .true. IF overlap is to be
! CRAY      LOGICAL variable set equal to .true. IF running on a Cray
! SCHEME    INTEGER variable: =1 => molecular problem (Gauss_Legendre qu
!                             =2 => coulomb problem (Gauss_Laguerre quad
! CSTEST    LOGICAL , usually false.  true tests basis for mega=1
! LMATCH    LOGICAL , is set equal to true to WRITE out information
!           for the FEM_to_ABM or DVR_to_ABM transformation.
!-----------------------------------------------------------------------
OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile,form='formatted',status='old')
READ(In_Unit,NML=options,IOSTAT=IERR)
IF(IERR==-1)THEN
      CALL NameList_Default(Out_Unit, 'Options')
   ELSEIF(IERR/=0)THEN
      WRITE(Out_Unit,*)'ERROR with Namelist Options'
      WRITE(Msg_Unit,*)'ERROR with Namelist Options'
      STOP 'ReadABM: Options'
   ENDIF
WRITE(Out_Unit,options)
CLOSE(Unit=In_Unit)
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
!    subs(i) ...... name of the routine to turn on debug WRITEs;
!                   must be a Hollerith CHARACTER string of no more than
!                   eight CHARACTERs
!    iprt(i) ...... what level of debug WRITEs:
!                    0=NONE
!                    1=minimum
!                    2=moderate
!                    3=maximum
!                    4=combination (must READ in iter(j,i), j=1 to 3)
!   iter(j,i) ...... how many times to loop thru routine i with the
!                    jth debug level turned on (j=1 to 3)
!-----------------------------------------------------------------------
OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile,form='formatted',status='old')
READ(In_Unit,NML=debug,IOSTAT=IERR)
IF(IERR==-1)THEN
      CALL NameList_Default(Out_Unit, 'Debug')
   ELSEIF(IERR/=0)THEN
      WRITE(Out_Unit,*)'ERROR with Namelist Debug'
      WRITE(Msg_Unit,*)'ERROR with Namelist Debug'
      STOP 'ReadABM: Debug'
   ENDIF
CLOSE(Unit=In_Unit)
CALL DebugOptions(Out_Unit) !WRITE(Out_Unit,NML=debug)
!------------------------------------------------------------------
! Determine printing options.
!------------------------------------------------------------------
CALL popt('ReadABM ', little, medium, full, ithcall, ithsub)
IF(little)THEN
  WRITE(Out_Unit,options)
  WRITE(Out_Unit,*)'$DEBUG'
  WRITE(Out_Unit,*)"In routine ReadABM"
  WRITE(Out_Unit,*)' NSUBS=',nsubs
  DO i = 1, nsubs
    WRITE(Out_Unit,'(" IPRT=",i4,3x,"SUBS=",A)')iprt(i), TRIM(subs(i))
    IF(iprt(i)==4)THEN
      WRITE(Out_Unit,*)' ITER=',iter(1,i),iter(2,i),iter(3,i)
    ENDIF
  ENDDO
  WRITE(Out_Unit,*)'$END'
ENDIF

!----------------------------------------------------------------------
! Specify the following:
!   JTOT    total angular momentum
!   PARITY  parity (0=even-parity, 1=odd-parity)
!   SYMMETRYis true only IF the diatom in that channel is homonuclear.
!           SYMMETRY can be true for only one arrangement channel.
!   JEVEN   true is even j states of the asymptotic diatoms
!           false is odd j states of the asymptotic diatoms
!----------------------------------------------------------------------
INQUIRE(File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile,Exist=There)
IF(There)THEN
OPEN(In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile, form='formatted',status='old')
READ(In_Unit,NML=momentum,IOSTAT=IERR)
IF(IERR==-1)THEN
      CALL NameList_Default(Out_Unit, 'Momentum')
   ELSEIF(IERR/=0)THEN
      WRITE(Out_Unit,*)'ERROR with Namelist Momentum'
      WRITE(Msg_Unit,*)'ERROR with Namelist Momentum'
      STOP 'ReadABM: Momentum'
   ENDIF
jtot = jref
megamax = jref
MxMega=MegaMax
WRITE(Out_Unit,NML=momentum)
symmetry = nsymc(1)
CLOSE(Unit=In_Unit)
ELSE
STOP "READABM"
ENDIF


!----------------------------------------------------------------------
! READ in the distances used to make the surface functions:
!
!      rhomin      Value of the first rho for calculating
!                  surface functions. (units of Bohr)
!      rhomax      Value of the last rho for calculating
!                  surface functions. (units of Bohr)
!      deltarho    fractional increase in rho at each step.
!----------------------------------------------------------------------
!OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile, form='formatted', status='old')
!READ(In_Unit,NML=sfundist,IOSTAT=IERR)
!IF(IERR==-1)THEN
!  CALL NameList_Default(Out_Unit, 'SFunDist')
!ELSEIF(IERR/=0)THEN
!  WRITE(Out_Unit,*)'ERROR with Namelist SFunDist'
!  WRITE(Msg_Unit,*)'ERROR with Namelist SFunDist'
!  STOP 'ReadABM: SFunDist'
!ENDIF
NRho=0
WRITE(Out_Unit,'("LastRho=",I5)')LastRho
DO IthRho=1,LastRho
   IF(Basis_Type(IthRho)=='ABM_APH')THEN
      IF(NRho==0)IthStart_ABM=IthRho
      NRho=NRho+1
      IthEnd_ABM=IthRho
   ENDIF
ENDDO
WRITE(Out_Unit,'("IthStart_ABM=",I5," IthEnd_ABM=",I5," NRho=",I5)')IthStart_ABM,IthEnd_ABM,NRho
numrho = nrho
IF(little) WRITE(Out_Unit,sfundist)
CLOSE(Unit=In_Unit)
WRITE(Out_Unit,NML=sfundist)

!  --------------------------------------------------------------------
!   eigmin  min eigenvalue of overlap matrix to be kept in making a
!           linearly independent basis.
!   ovrerr  tolerated difference in overlaps evaluated in different
!           arrangement channels.
!   ngood   number of accurate eigenvalues and eigenfunctions wanted.
!  --------------------------------------------------------------------
OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile, form='formatted',status='old')
   READ(In_Unit,NML=convrg,IOSTAT=IERR)
   !nsfunc=SUM(NFreq)
   IF(IERR==-1)THEN
      CALL NameList_Default(Out_Unit, 'Convrg')
   ELSEIF(IERR/=0)THEN
      WRITE(Out_Unit,*)'ERROR with Namelist Convrg'
      WRITE(Msg_Unit,*)'ERROR with Namelist Convrg'
      STOP 'ReadABM: Convrg'
   ENDIF
   IF(ngood<=1.or.ngood>nfreq(0))ngood=max(1,nfreq(0)/2)
   IF(nave1<=1)nave1=1
   IF(nave2<=1.or.nave2>ngood)nave2=max(1,ngood/2)
   WRITE(Out_Unit,NML=convrg)
CLOSE(Unit=In_Unit)


!----------------------------------------------------------------------
!   minvib    minimum vibrational quantum number for each arrangement
!             channel.
!   maxvib    maximum vibrational quantum number for each arrangement
!             channel.
!   jmin      minimum rotational quantum number for each vibrational
!             state and each arrangement, for mega=0,
!             jmin(nu,karran,mega).
!   jmax      maximum rotational quantum number for each vibrational
!             state and each arrangement channel.
!----------------------------------------------------------------------
OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile, form='formatted',status='old')
READ(In_Unit,NML=quantum,IOSTAT=IERR)
IF(IERR==-1)THEN
   CALL NameList_Default(Out_Unit, 'Quantum')
ELSEIF(IERR/=0)THEN
   WRITE(Out_Unit,*)'ERROR with Namelist Quantum'
   WRITE(Msg_Unit,*)'ERROR with Namelist Quantum'
   STOP 'ReadABM: QUantum'
ENDIF
CLOSE(In_Unit)
DO iarran=1,narran
   nlegndre(iarran)=jmax(0,iarran)
   DO Vib=MinVib(iarran),MaxVib(iarran)
      jmax(Vib,iarran)=jmax(0,iarran)
   ENDDO
ENDDO
WRITE(Out_Unit,NML=quantum)

ALLOCATE(MxBasisMega(0:MegaMax), MxBasissMega(0:MegaMax))
DO Mega=0,MegaMax
   MxBasisMega(Mega)=0
   DO iarran=1,narran
      DO vib=MinVib(iarran),MaxVib(iarran)
         MxBasisMega(Mega)=MxBasisMega(Mega)+(jmax(vib,iarran)+1)
      ENDDO
   ENDDO
   WRITE(Out_Unit,*)'For Mega=',Mega,'  MxBasisMega(Mega)=',MxBasisMega(Mega)
ENDDO
MxBasis =MAXVAL(MxBasisMega)
MxBasiss=MxBasis
!----------------------------------------------------------------------
! IF symmetry decoupling is desired make sure quantum numbers are correct.
! First check vibrational quantum numbers.
!----------------------------------------------------------------------
IF(symmetry)THEN
   IF(minvib(2)/=minvib(3).or.maxvib(2)/=maxvib(3))THEN
      WRITE(Out_Unit,*)'The number of vibrational ', 'states must be ', 'equivalent for channels with identical atoms'
      STOP 'ReadABM'
   ENDIF

!----------------------------------------------------------------------
! Now check rotational quantum numbers.
!----------------------------------------------------------------------
   DO  vib = minvib(2), maxvib(2)
      IF(jmin(vib,2, 0)/=jmin(vib, 3, 0).or.jmax(vib,2)/=jmax(vib,3))THEN
         WRITE(Out_Unit,*)'The number of rotational states for ', &
            'each vibrational state must be equal ', 'for the channels with identical atoms'
         STOP 'ReadABM'
      ENDIF
   ENDDO

   DO vib = minvib(1), maxvib(1)
      IF(jeven)THEN
         IF(mod(jmin(vib,1, 0),2)/=0 .or.mod(jmax(vib,1),2)/=0)THEN
            WRITE(Out_Unit,*)'jmin and jmax must be even for ',  'symmetry channel'
            STOP 'ReadABM'
         ENDIF
      ELSE
         IF(mod(jmin(vib,1, 0),2)/=1 .or.mod(jmax(vib,1),2)/=1)THEN
            WRITE(Out_Unit,*)'jmin and jmax must be odd for ', 'symmetry channel'
            STOP 'ReadABM'
         ENDIF
      ENDIF
   ENDDO

   DO Mega=0,MegaMax
      MxBasisMega(Mega)=0
      MxBasissMega(Mega)=0
      DO iarran=1,narran
         DO vib=MinVib(iarran),MaxVib(iarran)
             IF(iarran==1)THEN
               MxBasisMega(Mega) =MxBasisMega(Mega) +(jmax(vib,iarran)+1)
               MxBasissMega(Mega)=MxBasissMega(Mega)+(jmax(vib,iarran)/2+1)
            ELSEIF(iarran==3)THEN
               MxBasisMega(Mega)=MxBasisMega(Mega)+(jmax(vib,iarran)+1)
            ELSE
               MxBasisMega(Mega) =MxBasisMega(Mega) +(jmax(vib,iarran)+1)
               MxBasissMega(Mega)=MxBasissMega(Mega)+(jmax(vib,iarran)+1)
            ENDIF
         ENDDO
      ENDDO
      WRITE(Out_Unit,*)'For Mega=',Mega,'  MxBasisMega(Mega)=',MxBasisMega(Mega)
   ENDDO
   MxBasis =MAXVAL(MxBasisMega)
   MxBasiss=MAXVAL(MxBasissMega)
ENDIF
Mxl=MAXVAL(jmax)+1
Mxn=MAXVAL(maxvib)+1

IF(little) WRITE(Out_Unit,quantum)


!----------------------------------------------------------------------
!  NHERMT   number of Gauss_Hermite points used in
!           calculating the Hamiltonian matrix elements.
!  NGLEGN   number of Gauss_Legendre quadrature points.  nglegn(1) must
!           be even IF symmetry is true.
!  NOSCIL   no. of harmonic oscillator basis fcns in each arrangement.
!  NLEGNDRE number of legendre polynomials used in basis.
!  RE       equilibrium position of the diatom for each arrangement
!           channel in au.
!  RX       RX times RE is the position of the minimum  of the
!           parabola which defines the harmonic oscillator basis.
!           This parameter can usually be set equal to 1.1 for each
!           arrangement channel.  This parameter should be greater than
!           one because the long range part of the potential is softer
!           than short range part of the potential.
!  WEAU     fundamental frequency of diatomic in a.u.
!  CALPHA   asymptotic scaling factor for steepness of gaussians.
!           It should be near unity. Has effect of scaling weau.
!  WEXEAU   anharmonicity constant of diatomic in a.u.
!  ANHARM   asymptotic adjustment factor for anharmonicity.  Near unity,
!           empirical.  IF wexeau or anharm is 0, uses harmonic basis.
!  DALPHA   depth of Morse potential fit to scaling factor.
!  BALPHA   steepness parm of Morse potential fit to scaling factor.
!  RALPHA   equil parm of Morse potential for scaling factor.  Near
!           transition state rho.  These three params allow the
!           gaussians to spread in transition region and
!           shrink at small rho.
!  INTWT    weighting of each arrangement channel in doing integrals.
!           Usually 1 for each channel.  For very tight channels with
!           good quadratures, set it to 2.  THEN, all quadratures
!           involving that arrangement are done in that arrangement.
!-----------------------------------------------------------------------
INQUIRE(File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile,Exist=There)
IF(There)THEN
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
ELSE
   WRITE(Out_unit,*)"Required file does not exist, InputFile=", InputDIR(1:LEN(TRIM(InputDIR)))//InputFile
   STOP "ReadABM"
ENDIF

!  ---------------------------------------------------------------------
!  check nglegn for symmetric case
!  ---------------------------------------------------------------------
   IF(symmetry)THEN
      IF(mod(nglegn(1),2)/=0)THEN
         WRITE(Out_Unit,*)'nglegn(1) must be even ','for symmetric case'
         STOP 'ReadABM'
      ENDIF
   ENDIF
MxGlegn=MAXVAL(NGlegn)
MxHermt=MAXVAL(NHermt)
! Need to be deleted oliver
mxn=max(mxn,maxval(noscil))

IF(little) WRITE(Out_Unit,gauss)
nangch(1)=1
!-----------------------------------------------------------------------
! Determine MxAng
!-----------------------------------------------------------------------
DO iarran=1,narran
   IF(symmetry)THEN
      IF(iarran==1)THEN
         nangch(iarran+1) = nangch(iarran) + nhermt(iarran)*(nglegn(iarran)/2)
      ELSEIF(iarran==2)THEN
         nangch(iarran+1) = nangch(iarran) + nhermt(iarran)*nglegn(iarran) !GregParker  TEMPMOD ?
         !nangch(iarran+1) = nangch(iarran) + nhermt(iarran)*nglegn(iarran)/2
      ELSE
         nangch(iarran+1) = nangch(iarran) !GregParker  TEMPMOD ?
         !nangch(iarran+1) = nangch(iarran) + nhermt(iarran)*nglegn(iarran)/2
      ENDIF
   ELSE
      nangch(iarran+1) = nangch(iarran) + nhermt(iarran)*nglegn(iarran)
   ENDIF
ENDDO
WRITE(Out_Unit,*)'At pos1 nangch=',nangch
WRITE(Out_Unit,*)'At pos1 nhermt=',nhermt
WRITE(Out_Unit,*)'At pos1 nglegn=',nglegn


CLOSE(Unit=In_Unit)

OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile, status='old',form='formatted')
READ(In_Unit,NML=qcasep,IOSTAT=IERR)
IF(IERR==-1)THEN
      CALL NameList_Default(Out_Unit, 'QCasep')
   ELSEIF(IERR/=0)THEN
      WRITE(Out_Unit,*)'ERROR with Namelist QCasep'
      WRITE(Msg_Unit,*)'ERROR with Namelist QCasep'
      STOP 'ReadABM: QCasep'
   ENDIF
CLOSE(Unit=In_Unit)
WRITE(Out_Unit,NML=qcasep)
WRITE(Out_Unit,*)'qcase=',qcase
IF(qcase.AND.RhoStartABM_APH>qswitch)THEN
   WRITE(Out_Unit,*)'Changing qcase to false since RhoStartABM_APH>qswitch'
   WRITE(OUT_Unit,*)'RhoStartABM_APH=',RhoStartABM_APH,'  qswitch=',qswitch
   qcase=.false.
   WRITE(Out_Unit,*)'qcase=',qcase
ENDIF
IF(qcase)THEN
   thdiv(1) = 0.d0
   chdiv(1) = 0.d0
   OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile, status='old',form='formatted')
   READ(In_Unit,NML=qdivp,IOSTAT=IERR)
   IF(QCase)Chdiv(1)=-Pi/2
   IF(IERR==-1)THEN
      CALL NameList_Default(Out_Unit, 'QDivp')
   ELSEIF(IERR/=0)THEN
      WRITE(Out_Unit,*)'ERROR with Namelist QDivp'
      WRITE(Msg_Unit,*)'ERROR with Namelist QDivp'
      STOP 'ReadABM: QDivp'
   ENDIF
   WRITE(Out_Unit,NML=qdivp)
   CLOSE(Unit=In_Unit)
!         IF(thdiv(1)/=0.d0.or.chdiv(1)/=0.d0)THEN
!          WRITE(Out_Unit,*)'Error: thdiv(1)/=0 .and/or chdiv(1)/=0'
!          WRITE(Out_Unit,*)thdiv(1),chdiv(1)
!          STOP 'ReadABM'
!         ENDIF
   thdiv(nthdiv+1)=pi/2
   chdiv(1)=-pi/2
   chdiv(nchdiv+1)=pi/2
   DO idiv = 1, nchdiv+1
      chdiv(idiv) = chdiv(idiv)
   ENDDO
ENDIF

IF(qcase)THEN
   nth = 0
   DO ithdiv = 1, nthdiv
      nth = nth + nqth(ithdiv)
   ENDDO
   nch = 0
   DO ichdiv = 1, nchdiv
      nch = nch + nqch(ichdiv)
   ENDDO

   nang = nth*nch
   IF(symmetry)THEN
      IF(mod(nch,2)/=0)THEN
         WRITE(Out_Unit,*)'nch must be even for symmetry true'
         STOP 'setbasis'
      ENDIF
      nang=nang/2
   ENDIF
   nangch(1) = 1
   nangch(2) = 1 + nang
   nangch(3) = 1 + nang
   nangch(4) = 1 + nang
   WRITE(Out_Unit,*)'nang,nth,nch=',nang,nth,nch
ELSE
   nang = nangch(narran+1) - 1
ENDIF

MxAng=nang
!revisd by X.L. on 08-19-2006
maxtemp=0
DO i=1,narran
   IF(i==1.AND.symmetry)THEN
      maxtemp=maxtemp+nhermt(i)*nglegn(i)/2
   ELSE
      maxtemp=maxtemp+nhermt(i)*nglegn(i)
   ENDIF
ENDDO

mxang=max(nang,maxtemp)
! Print variable used to dimension arrays
WRITE(Out_Unit,*)'MxAng=   ', MxAng
WRITE(Out_Unit,*)'MxBasis= ', MxBasis
WRITE(Out_Unit,*)'MxBasiss=', MxBasiss
WRITE(Out_Unit,*)'MxGlegn= ', MxGlegn
WRITE(Out_Unit,*)'MxHermt= ', MxHermt
WRITE(Out_Unit,*)'Mxl=     ', Mxl
WRITE(Out_Unit,*)'Mxn=     ', Mxn
WRITE(Out_Unit,*)'MxMega   ', MxMega

DEALLOCATE(MxBasisMega, MxBasissMega)

RETURN
ENDSUBROUTINE ReadABM
