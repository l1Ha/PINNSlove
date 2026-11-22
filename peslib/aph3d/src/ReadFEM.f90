SUBROUTINE readfem(numrho)
USE Numeric_Kinds_Module
USE Numbers_Module
USE FileUnits_Module
USE PES_MODULE
!USE Converge_Module
USE Narran_Module
USE Parms_Module
USE nstate_Module
USE quad_Module
USE fbeta_Module
USE fast_Module
USE crays_Module
USE doenlvls_Module
USE opts_Module
USE SfnDis_Module
USE Boundary_Module
USE prntng_Module
USE Masses_Module
USE sublst_Module
USE elemnt_Module
USE scalag_Module
USE extra_Module
USE ngrid_Module
USE totj_Module
!!!USE gauss_Module !tempmod
USE pow_Module
USE mdata_Module
USE main_Module
USE Spectro_Module
USE Regins_Module
USE engpro_Module
USE Storage_Module
USE Integrat_Module
USE aphchl_Module
USE qlam_Module
USE propag_Module
USE Title1_Module
USE potzero_Module
USE rhowrt_Module
USE restart_Module
USE rhovalue_Module
USE csbasis_Module
USE Contour_Module
USE PlotNL_Module
USE Mshdat_Module
USE SavNam_Module
USE Stst_Module
USE InputFile_Module
!
!   $RCSfile: readfem.f,v $ $Revision: 1.15 $
!   $Date: 89/11/16 11:03:24 $
!   $State: Stable $
!

!            P U R P O S E   O F    S U B R O U T I N E

!-----------------------------------------------------------------------
! Routine to READ in all the input DATA
! written by b.j. archer
!-----------------------------------------------------------------------

IMPLICIT NONE
SAVE
LOGICAL iset  !, lsample
INTEGER i, ithcall, numrho, ithrhot, IErr
REAL(Kind=WP_Kind) deltarho, rhovals, ovrtol, deltamin, deltamax
EXTERNAL regons, popt, glegen, massfac
!------------------------------------------------------------------
! DATA statements
!------------------------------------------------------------------
DATA ithcall/0/
DATA deltamin/0.d0/
rho2nd = 9998.0
rhostop = 9999.0
!-----------------------------------------------------------------------
! Get the titles
!-----------------------------------------------------------------------
OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile, form='formatted',status='old')
   READ(In_Unit,NML=Title_Labels)
CLOSE(Unit=In_Unit)
WRITE(Out_Unit,'(A)')Title1
WRITE(Out_Unit,'(A)')Title2
WRITE(Out_Unit,'(A)')Title3
WRITE(Out_unit,*)
!--------------------------------------------------------------------
! READ in a label to uniquely identify the potential.
!--------------------------------------------------------------------
WRITE(Out_unit,'(A)') PES_Name
!-----------------------------------------------------------------------
!  Main option parameters:
! LMESHER   LOGICAL variable set equal to .true. IF mesher is to be
!           executed.
! LSFUNC    LOGICAL variable set equal to .true. IF sfunc is to be
!           executed.
! LMATELM   LOGICAL variable set equal to .true. IF matelm is to be
!           executed.
! LOVERLAP  LOGICAL variable set equal to .true. IF overlap is to be
!           executed.
! LPLOT     LOGICAL variable set equal to .true. for plotting
!           wavefunctions and potential.
! LAPH3D    LOGICAL variable set equal to .true. IF aph3d is to be
!           executed.
! LXSEC     LOGICAL variable set equal to .true. IF xsec is to be
!           executed.
! CRAY      LOGICAL variable set equal to .true. IF running on a Cray
! LENLVLS   LOGICAL variable set equal to .true. IF the propagator
!           is to calculate the energy correlations
! SCHEME    INTEGER variable: =1 => molecular problem (Gauss_Legendre qu
!                             =2 => coulomb problem (Gauss_Laguerre quad
!-----------------------------------------------------------------------
OPEN(Unit=In_Unit,File=TRIM(InputDIR)//InputFile, form='formatted',status='old')
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
      STOP 'ReadAPH_ph: Options'
   ENDIF
   WRITE(Out_Unit,NML=options)
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
OPEN(Unit=In_Unit,File=TRIM(InputDIR)//TRIM(InputFile),form='formatted',status='old')
   READ(In_Unit,NML=debug,IOSTAT=IERR)
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
   CALL DebugOptions(Out_Unit)   !WRITE(Out_Unit,NML=debug)
CLOSE(Unit=In_Unit)
!------------------------------------------------------------------
! Determine printing options.
!------------------------------------------------------------------
CALL popt('readfem ',little,medium,full,ithcall,ithsub)
IF(little)THEN
  WRITE(Out_unit,options)
  WRITE(Out_unit,*)'$DEBUG'
  WRITE(Out_Unit,*)"In routine ReadFEM"
  WRITE(Out_unit,*)' NSUBS=',nsubs
  DO 10 i=1,nsubs
    WRITE(Out_unit,16)iprt(i), TRIM(subs(i))
 16       FORMAT(" IPRT=",i4,3x,"SUBS=",A)
    IF(iprt(i)==4)THEN
      WRITE(Out_unit,*)' ITER=',iter(1,i),iter(2,i),iter(3,i)
    ENDIF
 10     CONTINUE
  WRITE(Out_unit,*)'$END'
ENDIF
!----------------------------------------------------------------------
! Specify the following:
!   JTOT   total angular momentum
!   IPAR   parity (IPAR=0=even-parity, IPAR=1=odd-parity)
!   NSYMC  is true only IF the diatom in that channel is homonuclear.
!          NSYMC can be true for only one arrangement channel.
!  NPAR   true means use parity decoupling for the 3-body system,
!         false means no parity decoupling
!  JEVEN  true is even parity states of the asymptotic diatoms
!         false is odd parity states of the asymptotic diatoms
!----------------------------------------------------------------------
OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile,form='formatted',status='old')
READ(In_Unit,NML=momentum,END=999)
jtot = jref
megamax = jref
WRITE(Out_unit,NML=momentum)
CLOSE(Unit=In_Unit)

!----------------------------------------------------------------------
! Set the zero of the potential (in Hartree's), IF needed.
!----------------------------------------------------------------------
OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile,form='formatted',status='old')
READ(In_Unit,  NML=potzero, IOSTAT=IERR)
IF(IERR==-1)THEN
   CALL NameList_Default(Out_Unit, 'PotZero')
   vzero=0.d0
ELSEIF(IERR/=0)THEN
   WRITE(Out_Unit,*)'ERROR with Namelist PotZero'
   WRITE(Msg_Unit,*)'ERROR with Namelist PotZero'
   STOP 'ReadAllData: PotZero'
ENDIF
WRITE(Out_Unit,NML=PotZero)
CLOSE(Unit=In_Unit)
!----------------------------------------------------------------------
! Save files to restart the code.
!
! ithsav    The files will be saved every ithsav rho value.  To turn off
!           this option set ithsav larger than the number of rho's wher
!           surface functions are calculated.
! lendir    Length of the name of the directory where the files
!           are to be saved on mass.
! ndir      A Hollerith string less than 70 CHARACTERs long
!           which contains the name of directory on mass where the
!           files are to be saved.
!----------------------------------------------------------------------
OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile, form='formatted',status='old')
READ(In_Unit,  NML=savnam, IOSTAT=IERR)
IF(IERR==-1)THEN
   CALL NameList_Default(Out_Unit, 'SavNam')
   ithsav=  100
ELSEIF(IERR/=0)THEN
   WRITE(Out_Unit,*)'ERROR with Namelist SavNam'
   WRITE(Msg_Unit,*)'ERROR with Namelist SavNam'
   STOP 'ReadAllData: SavNam'
ENDIF
WRITE(Out_Unit,NML=savnam)
IF(little)THEN
  WRITE(Out_unit,*)' $SAVNAM ITHSAV=',ithsav,', LENDIR=',lendir,','
  WRITE(Out_unit,410)ndir
 410    FORMAT(2x,'NDIR=',a70)
  WRITE(Out_unit,*)' $END'
ENDIF
CLOSE(Unit=In_Unit)
!----------------------------------------------------------------------
! READ in the distances used to make the surface functions:
!
!      rhomin      Value of the first rho for calculating
!                  surface functions. (units of Bohr)
!      rhomax      Value of the last rho for calculating
!                  surface functions. (units of Bohr)
!      deltarho    Length of the first sector. (unit of Bohr)
!      ithrho      IF ithrho is less than or equal to 1 a new
!                  calculation will be started.  IF ithrho is
!                  greater than 1 it is assumed that ithrho
!                  steps have already been calculated
!                  and the calculations will start at the ithrho+1
!                  step.
!       rhowrit    The surface functions are written out for use by
!                  the propagator after this rho value.
!       rhoend     IF restarting (ithrho>1), rhoend is the last
!                  good rho value.
!----------------------------------------------------------------------
!OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile, form='formatted',status='old')
!READ(In_Unit,  NML=sfundist, IOSTAT=IERR)
!IF(IERR==-1)THEN
!   CALL NameList_Default(Out_Unit, 'SFunDist')
!   nrho= 301
!   rhomin=    2.1d0
!   deltarho1= 0.05d-00        !5.478743543876281D-01
!   deltarho2= 0.00d-00        !1.0000000000000D-02,
!   rhomax=   8.10D0
!ELSEIF(IERR/=0)THEN
!   WRITE(Out_Unit,*)'ERROR with Namelist SFunDist'
!   WRITE(Msg_Unit,*)'ERROR with Namelist SFunDist'
!   STOP 'ReadAllData: SFunDist'
!ENDIF
WRITE(Out_Unit,NML=sfundist)
IF(little)THEN
  WRITE(Out_unit,*)' $SFUNDIST'
  WRITE(Out_unit,*)'  RHOMIN=',rhomin,' RHOMAX=',rhomax
  WRITE(Out_unit,*)'  DELTARHO1=',deltarho1,'   DELTARHO2=',deltarho2
  WRITE(Out_unit,*)'  NRHO=',nrho
  WRITE(Out_unit,*)' $END'
ENDIF
NRho=0
WRITE(Out_Unit,'("LastRho=",I5)')LastRho
DO IthRhot=1,LastRho
   IF(Basis_Type(IthRhot)=='FEM')THEN
      IF(NRho==0)IthStart_FEM=IthRhot
      NRho=NRho+1
      IthEnd_FEM=IthRhot
   ENDIF
ENDDO
LastRho=0
ithrho=1
nrho=IthEnd_FEM-IthStart_FEM+1
numrho = nrho
deltarho = deltarho1  !tempmod
WRITE(Out_unit,*) 'ithmax =',ithmax !tempmod
IF(rhostop > 0.d0 .AND. rhostop < 9997.d0)THEN !tempmod
  WRITE(Out_unit,*) 'rho2nd, rhostop =', rho2nd, rhostop !tempmod
ENDIF !tempmod
!IF(rhowrit>9998.d0)rhowrit=0.95d0*rhomax !tempmod
CLOSE(Unit=In_Unit)
!----------------------------------------------------------------------
! READ in the mesher DATA:
! NLT      initial number of elements in the theta direction.
! NLC      initial number of elements in the chi direction for each
!          arrangement channel.
! ENGEV    maximum total energy (units of eV) for scattering.
! NDIV     number of times each element can be subdivided.
! NFREQ    number of eigenenergies desired.
!  NXTRA   number of extra vectors used in sspace.
! NQ       number of basis vectors used (NQ=NXTRA+NFREQ).
! NODMAX   maximum number of nodal points desired.
! TOLEXP   exponent to determine feathering of the mesh IF scheme=2.
!  tolfn1   wavefunction error tolerance at first rho.
!  tolfn2   wavefunction error tolerance at subsequent rho values.
!  xkmesh  maximum wavenumber(au).  used in determining mesh.
!----------------------------------------------------------------------
IF(lmesher)THEN
   OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile,form='formatted',status='old')
   nxtra=0
   aterm=5.d+2
   bterm=1.d+5
   READ(In_Unit,NML=mshdat)
   WRITE(Out_Unit,NML=mshdat)
   IF(nave1_FEM<=0)nave1_FEM=1
   IF(nave2_FEM>nfreq(0))nave2_FEM=nfreq(0)
   IF(nave2_FEM<=0)nave2_FEM=nfreq(0)
   READ(In_Unit,NML=stst,END=999)
   CLOSE(Unit=In_Unit)
   WRITE(Out_Unit,NML=stst)
   CALL glegen(nquad,xg,wgt,-One,One)
!  IF nxtra is not set in mshdat, use Bathe's recommended value
! for the subspace method, no extra surface funcitons are
! used by lancoz.
  IF(nlancz==0)THEN
    IF(nxtra==0) nxtra=min(nfreq(0), 8)
    nq=nfreq(0)+nxtra
  ELSE
    nq=nfreq(0)+nxtra
  ENDIF
  WRITE(Out_unit,*)'nfreq,nxtra,nq=',nfreq,nxtra,nq
  IF(little) WRITE(Out_unit,mshdat)
  IF(little) WRITE(Out_unit,stst)
!----------------------------------------------------------------------
! calculate the number of rows (ntheta) and columns (nchi) in the virtua
! grid.
!-----------------------------------------------------------------------
  ntheta=nlt*(2**ndiv)*2+1
  nchi=0
  DO i=1,3
     nchi=nchi+nlc(i)*(2**ndiv)*2+1
     IF(i/=1)nchi=nchi-1
  ENDDO
  IF(nsymc(1))nchi=(nchi+1)/2
  nchi=2*nchi-1
  WRITE(Out_unit,*)'ntheta, nchi=',ntheta,nchi
ENDIF
IF(lplot)THEN
!-----------------------------------------------------------------------
!  READ in contour plotting parameters.
!  For plotting contours only (not needed otherwise)
!-----------------------------------------------------------------------
  OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile,status='old')
  READ(In_Unit,NML=contour,END=998)
  CLOSE(Unit=In_Unit)
  WRITE(Out_unit,NML=contour)
!-----------------------------------------------------------------------
!  READ in 3-D plotting parameters.
!  For plotting 3-dimensional graphs using grafic (not used otherwise)
!-----------------------------------------------------------------------
  OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile,status='old')
  READ(In_Unit,NML=plot,END=997)
  WRITE(Out_unit,NML=plot)
ENDIF
RETURN
 997  WRITE(Out_unit,*)'No plot DATA'
STOP 'readfem'
 998  WRITE(Out_unit,*)'No contour DATA'
STOP 'readfem'
 999  WRITE(Out_unit,*)'READ to END of aphdat, missing a NAMELIST'
STOP 'readfem'
!--------------------------***END-readfem***----------------------------
ENDSUBROUTINE ReadFEM
