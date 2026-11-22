SUBROUTINE readinbb
!            P U R P O S E   O F    S U B R O U T I N E
!   This routine reads in data from UNIT=In_Unit
!            O U T P U T    A R G U M E N T S
!  jtot     total angular momentum
!  parity   parity
!this routine is called by:
!         intbas
!this routine calls
!         popt
!-----------------------------------------------------------------------
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

!            I N C L U D E    S T A T E M E N T S
USE fileunits_Module
USE Narran_Module
USE PES_MODULE
USE gaussb_Module
USE totj_Module
USE quantb_Module
USE InputFile_Module
USE Title1_Module

IMPLICIT NONE
EXTERNAL popt
!-----------------------------------------------------------------------
!  Obtain printing options.
!-----------------------------------------------------------------------
LOGICAL little, medium, full
INTEGER ithcall, ithsub, IErr, iarran, vib
DATA ithcall /0/, ithsub/0/
DATA little /.true./, medium/.false./, full/.false./
CALL popt ('readinbb', little, medium, full, ithcall, ithsub)

!-----------------------------------------------------------------------
! Get the titles
!-----------------------------------------------------------------------
WRITE(Out_Unit,*)'begining of routine readinbb'

OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile,form='formatted',status='old')
   READ(In_Unit,NML=Title_Labels)
CLOSE(Unit=In_Unit,status='keep')
WRITE(Out_Unit,'(A)')Title1
WRITE(Out_Unit,'(A)')Title2
WRITE(Out_Unit,'(A)')Title3
WRITE(Out_Unit,*)

!----------------------------------------------------------------------
! Specify the following:
!   JTOT    total angular momentum
!   PARITY  parity (IPAR=0=even-parity, IPAR=1=odd-parity)
!   NSYMC   is true only if the diatom in that channel is homonuclear.
!           NSYMC can be true for only one arrangement channel.
!   NPAR    true means use parity decoupling for the 3-body system,
!           false means no parity decoupling
!   JEVEN   true is even parity states of the asymptotic diatoms
!           false is odd parity states of the asymptotic diatoms
!----------------------------------------------------------------------
OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile,form='formatted',status='old')
   READ(In_Unit,  NML=momentum, IOSTAT=IERR)
   IF(IERR/=0)THEN
      WRITE(Out_Unit,*)'ERROR with Namelist momentum'
      WRITE(Msg_Unit,*)'ERROR with Namelist momentum'
      STOP 'ReadAllData: Momentum'
   ENDIF
   WRITE(Out_Unit,NML=momentum)
symmetry = nsymc(1)
CLOSE(Unit=In_Unit,status='keep')
WRITE(Out_Unit,NML=momentum)


!----------------------------------------------------------------------
!   minvib    minimum vibrational quantum number for each arrangement
!             channel.
!   maxvib    maximum vibrational quantum number for each arrangement
!             channel.
!   jmin      minimum rotational quantum number for each arrangement
!             channel and each vibrational state.
!   jmax      maximum rotational quantum number for each arrangement
!             channel and each vibrational state.
!----------------------------------------------------------------------
OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile,form='formatted',status='old')
   READ(In_Unit,  NML=quantum, IOSTAT=IERR)
   IF(IERR==-1)THEN
      CALL NameList_Default(Out_Unit, 'Quantum')
      minvib=[0,0,0]
      maxvib=[4,4,4]
      jmin=0
      jmax(0,1:3)=(Maxvib+[1,1,1])*6 ! Use (maxvib+1)*(jmax)
   ELSEIF(IERR/=0)THEN
      WRITE(Out_Unit,*)'ERROR with Namelist Quantum'
      WRITE(Msg_Unit,*)'ERROR with Namelist Quantum'
      STOP 'ReadAllData: Quantum'
   ENDIF
   DO iarran=1,narran
      nlegndre(iarran)=jmax(0,iarran)
      DO Vib=MinVib(iarran),MaxVib(iarran)
         jmax(Vib,iarran)=jmax(0,iarran)
      ENDDO
   ENDDO
   WRITE(Out_Unit,NML=quantum)
   CLOSE(In_Unit)


!----------------------------------------------------------------------
!  NHERMT   number of Gauss_Hermite points used in
!           calculating the Hamiltonian matrix elements.
!  NGLEGN   number of Gauss_Legendre quadrature points.
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
!          It should be near unity. Has effect of scaling weau.
!  WEXEAU   anharmonicity constant of diatomic in a.u.
!  ANHARM   asymptotic adjustment factor for anharmonicity.  Near unity,
!           empirical.  if wexeau or anharm is 0, uses harmonic basis.
!  DALPHA   depth of Morse potential fit to scaling factor.
!  BALPHA   steepness parm of Morse potential fit to scaling factor.
!  RALPHA   equil parm of Morse potential for scaling factor.  Near
!           transition state rho.  These three params allow the
!           gaussians to spread in transition region and
!           shrink at small rho.
!  INTWT    weighting of each arrangement channel in doing integrals.
!           Usually 1 for each channel.
!           For very tight channels with good
!           quadratures, set it to 2.
!           Then, all quadratures involving that
!           arrangement are done in that arrangement.
!-----------------------------------------------------------------------
OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile,form='formatted',status='old')
   READ(In_Unit,  NML=gauss,  IOSTAT=IERR)
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
!?      zeta=  [1.000000D0, 1.000000D0, 1.000000D0]
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
      STOP 'ReadAllData: Gauss'
   ENDIF
   WRITE(Out_Unit,NML=gauss)
CLOSE(In_Unit,status='keep')
noscil=maxvib+1
CLOSE(Unit=In_Unit)

WRITE(Out_Unit,*)'END of routine readinbb'
RETURN
END
