SUBROUTINE DVR
USE Numeric_Kinds_Module
USE Numbers_Module
USE Convrsns_Module
USE PES_MODULE
USE DVR_Module
USE DVRMasses_Module
USE OneD_Module
USE TotJ_Module
USE Prvios_Module
USE DVR2_Module
USE Boundary_Module
USE FileUnits_Module
USE Parms_Module
USE System_Module
USE Dipole_Module
USE DVRNL_Module
USE InputFile_Module
USE TwoD_Module,h2d=>csurf1,urev=>csurf2,urec2=>ovlp
!     SUBROUTINE dvc2venw(dvc2vout,tape6=dvc2vout)
!
!     2-d symmetry-adapted dvr program for c2v symmetry
!
!     DVC2VENW is IDENTICAL to DVC2VEN, except that most of the
!     WRITE statements have been ELIMINATED
!
!     calculate surface functions, eigenvalues,
!     potential matrix elements
!     for a RANGE OF RHO VALUES, as well as the overlap matrices
!     between the surface functions at these rho values,
!     using Pack's hyperspherical coordinate
!     system (CPL,108(1984),333; JCP,87(1987),3888).
!
!     The rho range may be divided into SEVERAL sectors, EACH with
!     DIFFERENT basis parameters.
!
!     A legendre polynomial
!     dvr is used in the theta coordinate, and a symmetry adapted
!     basis is used in the chi coordinate.
!     c2v symmetry is used--this program calculates one
!     representation at a time as stored in variable sym.
!
!     system of units used in calculation is eV for energy,
!     au for distance and
!     amu for mass. 2-d eigenvalues are printed out in eV and in scaled
!     unitless values defined by (2*mu*rho^2/(hbar^2))*E
!
!
! Program DVR_dipole
!
! Modified from dvc2vnew to calculate upper state and lower
! state matrix elements and overlaps and the dipole matrix
! element between upper and lower surface functions at a
! particular rho value.
! This dipole matrix element is defined
!
!         <phi_t(upper)| mu(electric) |phi_t'(lower)>.
!
!  -------------------------------------------------------------------
! Note: At present, the upper and lower surfaces must have the
!       same dvr symmetry. To make this more general will require
!       substantial modification of the code.
!
!                               M. Braunstein; Winter 1990,91
!  --------------------------------------------------------------------
!
! Modified to include restart option
!
!                          -M. Braunstein, Fall 91
! Input:
!   Unit 5: dvc2vien
!    more details are given later
!    dipole - LOGICAL switch to DO either upper state or upper and lower
!             states and dipole matrix elements between them
!    system - variables for upper and lower states
!    sysu - unique variables for upper state
!    pesnamu - name of potential for upper state
!    sysl - unique variables for lower state
!    pesnaml - name of potential for lower state
!    pesnamd - name of electronic dipole moment function
!
!    restart - LOGICAL switch IF true will enable restart option
!    irhost - number of rho where you want to restart
!              irhost must be odd !!
!
! Output:
!   Unit 6: dvc2vout
!        7: Pmatrx_DVR:    potential matrix elements
!       OvrLap_DVR_Unit: Ovrlp_DVR :    overlaps
!       23: sfelevlu  :    eigenvalues of upper surface functions
!       26: wavef_DVR :    wave function saved for upper state
!
! Notes:
!  the symbol "CC" will precede sections which have been added
!  to dvc2vnew to produce DVR_dipole
!  Intermediate output:
!     units 12,13 - surfun1,2; upper state surface functions(dvr)
!     units 15,16 - tthch1,2; transformation matrices for
!                   upper state ( restart info - never tested)
!     unit 14 - nevtot
!     unit 22 - sfu01a - fbr surface functions asymptotically
!
IMPLICIT NONE
INTEGER NTheta, NChi, Npma, Novr, Nrho, ithrho, i, nch, ii, ipma, Ierr
INTEGER iovr, isfn, nunita, nunitb, nunitta, nunittb, nunitc, nunitd, nrpt
INTEGER nevtot, itheta, krho, irho, ibasis, nctem, ichi, j, isec
INTEGER jotu, jotl, nseg, nsec, nrhos1, nevsum, nthnch, nreth, nchith
INTEGER m, n, l, ksf, ip, nchip, nthetap, jj, kk, nrhdum, nsfdum
INTEGER nthdum, nchdum, k, nsym, index, NTheta_Plt, NChi_Plt
REAL(KIND=WP_Kind) tdum, etime, tarray(2), tthen, deltat
REAL(Kind=WP_Kind) delrho, delpevm, asegu, bsegu, rhoend, rhofinal, rhodum
REAL(Kind=WP_Kind) rho, rhoz, rholast, rhonext, rhos, aseg, bseg
REAL(Kind=WP_Kind) rz, chii, eigenav, eigenav1, eigenav2, potmax
REAL(Kind=WP_Kind) pev2dmax, rhos1, dz, chi1z, dist, vrhoz, r, d, chi1
REAL(Kind=WP_Kind) vrhos, udvr, ucmat, ctmat, theta, ctuc, pevms1, CSURF1
REAL(Kind=WP_Kind) sumc, tchit, tthetat, CSURFBR, stuff, ASEGL, BSEGL
REAL(Kind=WP_Kind), ALLOCATABLE:: ThetaVal(:), ChiVals(:), Phi_Plt(:)
INTEGER dtvalues(8)
LOGICAL psurfonl, mid_zero, once
CHARACTER(LEN=6) SFType
CHARACTER(LEN=3) sym
CHARACTER(LEN=3) syu,syl
CHARACTER(LEN=10) hour
CHARACTER(LEN=8) today
CHARACTER(LEN=5) curzone
EXTERNAL anchica1,anchicb1,anchica2,anchicb2
EXTERNAL td2fbca1,td2fbcb1,td2fbca2,td2fbcb2
EXTERNAL anchicp,anchicpp,td2fbcp,td2fbcpp
!  -------------------------------------------------------------------
!  commons and their explanations follow
!  -------------------------------------------------------------------
!
!     the following contains much information about the
!     mass-scaled jacobi coordinates used to construct the
!     hyperspherical system. for this program, the only quantities
!     used are mtot,mass,scale and twomu. NOTE: the mass-scaling factor
!     c_alpha is defined here to be the inverse of that used by Pack
!     and Parker. This variables are determined in SUBROUTINE setup.
!
COMMON /waste/ stuff(nt2max*(nt2max+3)/2)
!
!     h2d--2-d hamiltonian
!     later, h2d is no longer needed, the final 2-d eigenvectors are stored in h2d
!
!COMMON /twod/ h2d(nt2max,nt2max),        urev(nt2max,nt2max),      urec2(nt2max,nt2max)
!
!  -------------------------------------------------------------------
!  DATA statements follow.
!  some of these DATA statements set default values that are only used
!  IF the DATA READ in is incomplete.
!  ------------------------------------------------------------------

DATA ntheta,nchi/10,10/,npma/1000/,novr/2000/
DATA sym/'a1 '/
!
!  -----------------------------------------------------------------
!  DIMENSION statements
!  -----------------------------------------------------------------
DIMENSION csurf1(nthetmax*nchimax*nsfnmax,1)
DIMENSION csurfbr(nsfnmax)
!
!     DIMENSION rhos(3),udvr(nthetmax,nchimax),urev(nt2max,nt2max)
DIMENSION rhos(3),udvr(nthetmax,nchimax)
DIMENSION ucmat(nchimax,nchimax),ctmat(nchimax,nchimax)
!     DIMENSION ctuc(nchimax,nchimax),urec2(nt2max,nt2max)
DIMENSION ctuc(nchimax,nchimax)
DIMENSION dist(3)
!
DIMENSION rhos1(15),pevms1(15)
DIMENSION aseg(15),bseg(15)
DIMENSION tchit(nchimax,nchimax),tthetat(nthetmax,nthetmax)
!
DIMENSION asegu(15), asegl(15), bsegu(15), bsegl(15)

!
!  -------------------------------------------------------------------
!  OPEN some of the input and output files
!  -------------------------------------------------------------------
WRITE(Out_Unit,*)'Called    DVR'
CLOSE(Out_Unit)
OPEN(Unit=Out_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'Output/DVR.txt',Form='FORMATTED')
!OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile,     form='formatted',   status='old')
OPEN(Unit=dvr11,File=OutDIR(1:LEN(TRIM(OutDIR)))//'Output/dvrpts.txt',form='formatted',   status='unknown')
OPEN(Unit=dvr14,File=OutDIR(1:LEN(TRIM(OutDIR)))//'Output/nevtot.txt',form='formatted',   status='unknown')
OPEN(Unit=dvr23,File=OutDIR(1:LEN(TRIM(OutDIR)))//'GraphicsOut/Sfelevl_DVR.rbw',form='formatted',   status='unknown')

OPEN(Unit=dvr12,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/Surfun1_DVR.bin',form='unformatted', status='unknown')
OPEN(Unit=dvr13,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/Surfun2_DVR.bin',form='unformatted', status='unknown')
OPEN(Unit=dvr15,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/Tthch1_DVR.bin', form='unformatted', status='unknown')
OPEN(Unit=dvr16,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/Tthch2_DVR.bin', form='unformatted', status='unknown')
OPEN(Unit=dvr22,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/Sfunfbr_DVR.bin',form='unformatted', status='unknown')
!
!  -------------------------------------------------------------------
!     begin executable statements
!  -------------------------------------------------------------------
mid_zero = .true.
!  -------------------------------------------------------------------
!   Determine the time and today's date and WRITE them to unit Out_Unit.
!  ------------------------------------------------------------------

CALL Date_And_Time(today, hour, curzone, dtvalues)
CALL Print_Date_And_Time(Today, Hour, CurZone, DtValues, Out_Unit)
!  ------------------------------------------------------------------
!   WRITE headers
!  -----------------------------------------------------------------

WRITE(Out_Unit,'(/" A general 2-d c2v dvr code , june 1988"/)')
WRITE(Out_Unit,'("with modifications to calculate dipole matrix elements added 1991"/)')

!  -------------------------------------------------------------------
!     READ in input
!  -------------------------------------------------------------------
OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile,     form='formatted',   status='old')
   READ(In_Unit, NML=rstart, IOSTAT=IERR)
   WRITE(Out_Unit,NML=rstart)
CLOSE(In_Unit)

OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile,     form='formatted',   status='old')
   READ(In_Unit, NML=dipole, IOSTAT=IERR)
   WRITE(Out_Unit,NML=dipole)
CLOSE(In_Unit)

OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile,     form='formatted',   status='old')
   READ(In_Unit, NML=system)    !TMPModGregParker
   WRITE(Out_Unit,NML=system)    !TMPModGregParker
CLOSE(In_Unit)

OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile,     form='formatted',   status='old')
   READ(In_Unit, NML=sysu, IOSTAT=IERR)
   WRITE(Out_Unit,NML=sysu)
CLOSE(In_Unit)

IF(nsurf>1)THEN
   READ(In_Unit, NML=sysl, IOSTAT=IERR)
   WRITE(Out_Unit,NML=sysl)
ELSE
   PesNamu=PES_Name
ENDIF

NRho=0
WRITE(Out_Unit,'("LastRho=",I5)')LastRho
DO IthRho=1,LastRho
   IF(Basis_Type(IthRho)=='DVR')THEN
      IF(NRho==0)IthStart_DVR=IthRho
      NRho=NRho+1
      IthEnd_DVR=IthRho
   ENDIF
ENDDO
NRhoMx=NRho
WRITE(Out_Unit,'("IthStart_DVR=",I5," IthEnd_DVR=",I5," NRho=",I5)')IthStart_DVR,IthEnd_DVR,NRho

   11 FORMAT(a25)
!
!  -------------------------------------------------------------------
! WRITE out input to unit 6
!  -------------------------------------------------------------------
!
WRITE(Out_Unit,rstart)
WRITE(Out_Unit,'(/"input echo data for both states"/)')
WRITE(Out_Unit,'(A,3(ES22.14))') ' mass = ',(mass(i),i=1,3)
WRITE(Out_Unit,'(A,3(I22))') ' nthesec = ',(nthesec(i),i=1,Num_DVR_Sectors)
WRITE(Out_Unit,'(A,3(I22))') ' nchisec = ',(nchisec(i),i=1,Num_DVR_Sectors)
WRITE(Out_Unit,'(A,3(ES22.14))') ' rendsec = ',(rendsec(i),i=1,Num_DVR_Sectors)
WRITE(Out_Unit,*)  'peigv = ',peigv
WRITE(Out_Unit,*)  'debug = ',debug
WRITE(Out_Unit,*)  'medium = ',medium
WRITE(Out_Unit,'(A,3(ES22.14))') ' potmsec = ',(potmsec(i),i=1,Num_DVR_Sectors)
WRITE(Out_Unit,'(A,3(ES22.14))') ' pevmsec = ',(pevmsec(i),i=1,Num_DVR_Sectors)
WRITE(Out_Unit,'(A,3(I22))') ' ltheta = ',ltheta
WRITE(Out_Unit,'(A,3(ES22.14))') ' rho1 = ',rho1
WRITE(Out_Unit,*)  'matelem = ',matelem
WRITE(Out_Unit,*)  'h3sys = ',h3sys
WRITE(Out_Unit,'(A,3(ES22.14))') ' rhofbr1 = ',rhofbr1
WRITE(Out_Unit,'(A,3(ES22.14))') ' rhofbr2 = ',rhofbr2
WRITE(Out_Unit,'(A,3(I22))') ' nfbr = ',nfbr
WRITE(Out_Unit,'(A,3(I22))') ' nsfunc = ',nsfunc
WRITE(Out_Unit,'(A,3(I22))') ' nrow = ',nrow
WRITE(Out_Unit,'(A,3(I22))') ' ncol = ',ncol
WRITE(Out_Unit,'(A,3(I22))') ' ngood1 = ',ngood1
WRITE(Out_Unit,'(A,3(I22))') ' ngood2 = ',ngood2
WRITE(Out_Unit,'(A,3(I22))') ' Num_DVR_Sectors = ',Num_DVR_Sectors
WRITE(Out_Unit,'(A,3(I22))') ' initsec = ',initsec
WRITE(Out_Unit,'(A,3(I22))') ' nrhomx = ',nrhomx
WRITE(Out_Unit,'(A,3(ES22.14))') ' fmpotma = ',fmpotma
WRITE(Out_Unit,'(A,3(ES22.14))') ' deltar = ',deltar
WRITE(Out_Unit,*)  ' logspace = ',logspace
!
WRITE(Out_Unit,'(/"upper state input echo"/)')
WRITE(Out_Unit,'(/1x,a25/)')pesnamu
WRITE(Out_Unit,*)  'jevenu = ',jevenu
WRITE(Out_Unit,'(A,3(I22))') ' nrhos1u = ', nrhos1u
WRITE(Out_Unit,'(A,3(ES22.14))') ' rhos1u = ', (rhos1u(i),i=1,nrhos1u)
WRITE(Out_Unit,'(A,3(ES22.14))') ' pevms1u = ',(pevms1u(i),i=1,nrhos1u)
WRITE(Out_Unit,'(A,3(I22))') ' iparu = ',iparu
WRITE(Out_Unit,'(A,3(I22))') ' jtotu= ',jtotu
WRITE(Out_Unit,'(A,3(I22))') ' megau= ',megau
WRITE(Out_Unit,*) ' cstestu = ',cstestu
WRITE(Out_Unit,*)  ' pointu = ',pointu

!
IF(nsurf>1)THEN
   WRITE(Out_Unit,'(/"lower state input echo"/)')
   WRITE(Out_Unit,'(/1x,a25/)')pesnaml
   WRITE(Out_Unit,*)  'jevenl = ',jevenl
   WRITE(Out_Unit,'(A,3(I22))') ' nrhos1l = ', nrhos1l
   WRITE(Out_Unit,'(A,3(ES22.14))') ' rhos1l = ', (rhos1l(i),i=1,nrhos1l)
   WRITE(Out_Unit,'(A,3(ES22.14))') ' pevms1l = ',(pevms1l(i),i=1,nrhos1l)
   WRITE(Out_Unit,'(A,3(I22))') ' iparl = ',iparl
   WRITE(Out_Unit,'(A,3(I22))') ' jtotl= ',jtotl
   WRITE(Out_Unit,'(A,3(I22))') ' megal= ',megal
   WRITE(Out_Unit,*) ' cstestl = ',cstestl
   WRITE(Out_Unit,*)  ' pointl = ',pointl
   WRITE(Out_Unit,*) ' rhoendl = ',rhoendl
   WRITE(Out_Unit,'(/1x,a25/)') pesnamd
ENDIF


CLOSE(Unit=In_Unit)
!
 94   FORMAT(2x,"nthetmax=",i4,2x,"nchimax= ",i4,2x,"nt2max= ",i4)
 96   FORMAT(2x,"fmpotma=",f8.4/)
 97   FORMAT(2x,"Num_DVR_Sectors=",i2,2x,"initsec=",i2,2x,"nrhomx=",i4)
 98   FORMAT(2x,"nsfunc=",i4,2x,"nrow=",i3,2x,"ncol=",i3)
!
!  -------------------------------------------------------------------
!     set up mass-scaled coordinates. nch must always be set to
!     1 for some of the coordinate transformations to work.
!     66.9898927 is the value of 1/hbar^2 in the eV for energy-
!     au for distance-amu for mass system of units.
!     It equals 1822.888506/autoev (1986 values).
!  -------------------------------------------------------------------
!
nch=1
CALL setup2(.false.,one/autoev)   !3.67493d-2)
ndvptmax=nsfnmax            !TEMPMod
WRITE(Out_Unit,*)"Allocating Arrays: ndvptmax,nsfnmax=",ndvptmax,nsfnmax
ALLOCATE(h2d(ndvptmax,nsfnmax), urev(ndvptmax,nsfnmax), urec2(nsfnmax,nsfnmax))
!
IF( TRIM(pointg)=='cs' .AND. (mod(nchisec(1),2))/=0 )THEN
   WRITE(Out_Unit,*) 'nchisec(1) = ',nchisec(1)
   WRITE(Out_Unit,*) 'nchi must be a multiple of 2 for cs symmetry'
   STOP
ENDIF
!
DO ii=1,NSurf
   IF(ii==1)THEN
      pointg=pointu
      parity=iparu
      jtot=jtotu
      mega=megau
      cstest=cstestu
      jeven=jevenu
   ELSE
      pointg=pointl
      parity=iparl
      jtot=jtotl
      mega=megal
      cstest=cstestl
      jeven=jevenl
   ENDIF
   !  -----------------------------------------------------------------
   !   determine irreducible representation of functions wanted
   !  -----------------------------------------------------------------
   CALL repset(pointg, jeven, mega, parity, nsym)
   IF(ii==1)THEN
      IF(nsym==1) syu='a1 '
      IF(nsym==2) syu='a2 '
      IF(nsym==3) syu='b1 '
      IF(nsym==4) syu='b2 '
      IF(nsym==5) syu='ap '
      IF(nsym==6) syu='app'
   ELSE
      IF(nsym==1) syl='a1 '
      IF(nsym==2) syl='a2 '
      IF(nsym==3) syl='b1 '
      IF(nsym==4) syl='b2 '
      IF(nsym==5) syl='ap '
      IF(nsym==6) syl='app'
   ENDIF
   !  -------------------------------------------------------------------
   !  determine parameters of linear fit to pev2max
   !  on each segment where they differ.
   !  -------------------------------------------------------------------
   IF(ii==1)THEN
      DO nseg=1,nrhos1u-1
         delpevm=pevms1u(nseg)-pevms1u(nseg+1)
         delrho=rhos1u(nseg)-rhos1u(nseg+1)
         asegu(nseg)=delpevm/delrho
         bsegu(nseg)=pevms1u(nseg)-asegu(nseg)*rhos1u(nseg)
      ENDDO
   ELSE
      DO nseg=1,nrhos1l-1
         delpevm=pevms1l(nseg)-pevms1l(nseg+1)
         delrho=rhos1l(nseg)-rhos1l(nseg+1)
         asegl(nseg)=delpevm/delrho
         bsegl(nseg)=pevms1l(nseg)-asegl(nseg)*rhos1l(nseg)
      ENDDO
   ENDIF
ENDDO
!
!
! FUNCTIONAL VERSION I
! in the FIRST rho segment (nsec=1) only, PEV2DMAX and POTMAX are
! chosen to be EXPONENTIAL FUNCTIONS of RHO: ap*exp(-alfap*rho).
! ap and alfap are calculated now.
!
!     alfap=(1.0d0/(rendsec(1)-rho1))*alog(pev2dma1/pev2dma2)
!     ap=pev2dma1*exp(alfap*rho1)
!     WRITE(Out_Unit,95) alfap,ap
! 95   FORMAT(2x,'alfap=',f10.5,2x,'ap=',f10.5/)
!
! FUNCTIONAL VERSION II
! in the FIRST rho segment (nsec=1) only, PEV2DMAX and POTMAX are
! chosen to be HYPERBOLIC FUNCTIONS of RHO: ap/(rho**alfap)
! ap and alfap are calculated now
!
!     alfap=dlog(pev2dma1/pev2dma2)/dlog(rendsec(1)/rho1)
!     ap=pev2dma1*(rho1**alfap)
!     WRITE(Out_Unit,95) alfap,ap
!
nsec=initsec-1
isec=0
psurfonl=.false.
rhoend=rendsec(initsec)
rhofinal=rendsec(Num_DVR_Sectors)
!
!
! Determine the number of rho values to be stored in a given
! pmaXX and ovrXX file, respectively.
!
!  Activate these calculations for multiple pmatrix/Ovrlp files
!     npma = 420000/(2*nsfunc*nsfunc)
!     novr = 2*npma
WRITE(Out_Unit,'(2x,"npmatrx=",i7,2x,"novrlp=",i7)') npma,novr
ipma = 1
iovr = 1
isfn = 1
!
IF(.NOT.restart) irhost=1
!
!  *********************************************************
! BEGIN MAIN LOOP OVER RHO VALUES AT WHICH SURFACE EIGENVALUES
! AND EIGENFUNCTIONS, POTENTIAL COUPLING AND OVERLAP MATRICES,
! ARE CALCULATED
!  ************************************************************
!
IF(MOD(nbreak,2)/=0)THEN
   WRITE(Msg_Unit,*)'nbreak must be even and greater than 0', nbreak
   WRITE(Out_Unit,*)'nbreak must be even and greater than 0', nbreak
   STOP 'ERROR: Stopping in DVR'
ENDIF
DO nrho=irhost,nrhomx
   !
   ! OPEN appropriate units
   !
   IF(nrho==irhost)THEN
      ! matrix elements
      OPEN(Unit=dvr7,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/PotMatrx_DVR.bin',form='unformatted',status='unknown')
      !OPEN(Unit=dvr7,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/'//pmanu(icoun),form='unformatted',status='unknown')
      !OPEN(Unit=dvr8,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/'//pmanl(icoun),form='unformatted',status='unknown')
      ! overlaps
      OPEN(Unit=dvr18,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/Ovrlp_DVR.bin',form='unformatted',status='unknown')
      !OPEN(Unit=dvr19,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/'//ovrnu(icoun),form='unformatted',status='unknown')
      !OPEN(Unit=dvr20,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/'//ovrnl(icoun),form='unformatted',status='unknown')
      ! dipole matrix elements
      !OPEN(Unit=dvr25,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/'//tmutp(icoun),form='unformatted',status='unknown')
      ! wavefunctions
      !OPEN(Unit=dvr26,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/'//wavu(icoun),form='unformatted',status='unknown')
      !OPEN(Unit=dvr27,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/'//wavl(icoun),form='unformatted',status='unknown')
   ENDIF
   !
   ! select rho value and shuffle nunits which contain the overlaps
   !
   IF(nrho==irhost)THEN
      rho=rho1
      nunita=dvr12
      nunitb=dvr13
      nunitta=dvr15
      nunittb=dvr16
      nunitc=dvr32
      nunitd=dvr33
   ELSE
      IF(logspace)THEN
         rho=(1.0d0 + deltar)*rho
      ELSE
         rho=rho + deltar
      ENDIF
   ENDIF
   WRITE(Msg_Unit,*)
   WRITE(Msg_unit,'("CurrentSector=",I5," at Rho=",ES22.15,A25)')CurrentSector,Basis_Dist(CurrentSector),Basis_Type(CurrentSector)
   WRITE(Out_Unit,*)
   WRITE(Out_Unit,'("CurrentSector=",I5," at Rho=",ES22.15,A25)')CurrentSector,Basis_Dist(CurrentSector),Basis_Type(CurrentSector)
   Rho=Basis_Dist(CurrentSector)
   rhoz=rho
   IF(logspace)THEN
      rholast=rhoz/(1.0d0 + deltar)
      rhonext=(1.0d0 + deltar)*rhoz
   ELSE
      rholast=rhoz - deltar
      rhonext=rhoz + deltar
   ENDIF
   IF(CurrentSector/=1)THEN
      rholast=Basis_Dist(CurrentSector-1)
   ELSE
      rholast=Sector_Boundaries(CurrentSector)
   ENDIF
   !
   ! select three points within rho sector to use in evaluating matrix element
   !
   rhos(1)=Sector_Boundaries(CurrentSector)
   rhos(2)=rhoz
   rhos(3)=Sector_Boundaries(CurrentSector+1)
   
   CurrentSector=CurrentSector+1
   IF(rho>rhofinal)THEN
      rho=rhofinal
   ENDIF
   IF(nunita==dvr12.AND.nrho/=irhost)THEN
      nunita=dvr13
      nunitb=dvr12
      nunitta=dvr16
      nunittb=dvr15
      nunitc=dvr33
      nunitd=dvr32
   ELSE
      nunita=dvr12
      nunitb=dvr13
      nunitta=dvr15
      nunittb=dvr16
      nunitc=dvr32
      nunitd=dvr33
   ENDIF
   !  --------------------------------------------------------------------
   ! inner loop over upper and lower PE surfaces
   !  --------------------------------------------------------------------
   !ERROR BELOW jtotu and jtotl needs to be defined
   jotu=0
   jotl=0
   DO 1960 ii=1,nsurf
      IF(ii==1)THEN
         PES_Name=pesnamu
         parity=iparu
         jtot=jotu
         mega=megau
         cstest=cstestu
         pointg=pointu
         sym=syu
      ELSE 
         PES_Name=pesnaml
         parity=iparl
         jtot=jotl
         mega=megal
         cstest=cstestl
         pointg=pointl
         sym=syl
      ENDIF
      IF(mod(nrho,npma) == 0)THEN
         CLOSE(dvr7)
         ipma = ipma + 1
         !
         !  PAUSE 'mulitple ipmas'
         !
         !  This is for unformatted WRITEs.
         !        OPEN(Unit=dvr7,File=OutDIR(1:LEN(TRIM(OutDIR)))//pman(ipma),status='unknown',form='unformatted')
         !  This is for formatted WRITEs.
         !       OPEN(Unit=dvr7,File=OutDIR(1:LEN(TRIM(OutDIR)))//pman(ipma),status='unknown')
      ENDIF
      IF(mod(nrho,novr) == 0)THEN
         CLOSE(dvr18)
         iovr = iovr + 1
      ENDIF
      tdum = 0.d0
      tthen = tarray(1)
      IF(ii==1)THEN
         WRITE(Out_Unit,'(/,"****upper surface****",/,2x,"*** nrho =",i4,2x,"rho =",f17.13," ***")')nrho,rho
      ENDIF
!      WRITE(Out_Unit,1006) nrho,nunita,nunitb
!      WRITE(Out_Unit,1007) nrho,nunitta,nunittb
! 1006 FORMAT(2x,'nrho=',i4,2x,'nunita =',i3,2x,'nunitb =',i3/)
! 1007 FORMAT(2x,'nrho=',i4,2x,'nunitta=',i3,2x,'nunittb=',i3/)
!
!
! IF we are doing the lower surface and rho > rhoendl
! DO not calculate anything (GOTO END of loop)
!
IF(ii==2.AND.rho>rhoendl) GOTO 1960
!
!
!  -----------------------------------------------------------
! beginning of the loop to determine the parameters of
! the dvr and eigenvector basis in a given rho segment
! and set up dvr transformations for theta and chi
!  ---------------------------------------------------------
!
!
! took out IF block to DO more than one symmetry for multiple
! rhos
!      IF(nrho==1.or.rho>rhoend)THEN
! added IF block to DO upper and lower surface on first rho
! Assume that upper and lower surfaces have the same symmetry
!Need to modified to use multiple nsec in DVR code. By X.L. 12-11-2006
!IF(ii/=2.AND.nrho==irhost)THEN
IF(ii/=2.AND.(nrho==irhost.or.rho>rhoend))THEN
!
588 CONTINUE
nsec=nsec+1
rhoend=rendsec(nsec)
ntheta=nthesec(nsec)
nchi=nchisec(nsec)
IF(Ntheta>NthetMax.or.Nchi>NchiMax)THEN
   WRITE(*,*)"Ntheta or Nchi exceeds maximum allowed values"
   WRITE(*,*)"Ntheta=",Ntheta,"  NthetMax=",NthetMax
   WRITE(*,*)"Nchi=",Nchi,"  NchiMax=",NchiMax
   STOP "DVR"
ENDIF
IF(rho>rhoend)GOTO 588
!
  IF(nrho == 1)THEN
    IF(nsec /= 1)THEN
      IF( (rho<rendsec(nsec-1)) .or. (rho > rhoend) )THEN
         WRITE(Out_Unit,*) 'for nrho == 1: IF nsec /= 1, (nsec =', nsec,'),'
         WRITE(Out_Unit,*) 'THEN rendsec(nsec-1) < rho < rendsec(nsec)'
         WRITE(Out_Unit,*) 'rendsec(nsec-1),rho,rendsec(nsec) =',rendsec(nsec-1),rho,rendsec(nsec)
         STOP 'Error: Stopping in DVR'
      ENDIF
    ENDIF
  ENDIF
!
  IF( TRIM(pointg)=='cs' .AND. (mod(nchi,2))/=0 )THEN
     WRITE(Out_Unit,*) 'nchi = ',nchi
     WRITE(Out_Unit,*) 'nchi must be a multiple of 2 for cs symmetry'
     STOP
  ENDIF
!
potmax=potmsec(nsec)
pev2dmax=pevmsec(nsec)
!
isec=1
!
!     set up dvr transformations for chi
!
!
! For this version, where the two surfaces have the same symmetry,
! lets take the upper surface symmetries to determine the
! basis
  parity=iparu
  jtot=jtotu
  mega=megau
  cstest=cstestu
  pointg=pointu
  jeven=jevenu
  sym=syu
  PES_Name=pesnamu
!        parity=iparl
!        jtot=jtotl
!        mega=megal
!        cstest=cstestl
!        pointg=pointl
!        jeven=jevenl
!        sym=syl
!        PES_Name=pesnaml
!
IF(sym=='a1 ')THEN
   CALL symset(anchica1,td2fbca1,nchi,debug,sym)
ELSEIF(sym=='a2 ')THEN
   CALL symset(anchica2,td2fbca2,nchi,debug,sym)
ELSEIF(sym=='b1 ')THEN
   CALL symset(anchicb1,td2fbcb1,nchi,debug,sym)
ELSEIF(sym=='b2 ')THEN
   CALL symset(anchicb2,td2fbcb2,nchi,debug,sym)
ELSEIF(sym=='ap ')THEN
   CALL symset(anchicp,td2fbcp,nchi,debug,sym)
ELSEIF(sym=='app')THEN
   CALL symset(anchicp,td2fbcpp,nchi,debug,sym)
ENDIF
!
!     set up dvr transformation for theta
!
CALL tdv2fb1(ttheta,pttheta,whtheta,ntheta,ltheta, nthetmax,debug)
!
!     note that number of basis functions is ntheta-ltheta+1 for theta
!     Therefore, we readjust the value of
!     ntheta before we start the REAL calculation.
!
!
! DO this once and once only
!
  ntheta=ntheta-ltheta+1
!
!
!     angmom determines the dvr j^2 matrix for theta and
!     puts it in antothet
!
CALL anmomthe(ttheta,pttheta,ntheta,nthetmax,ltheta,antothet,thetaj, debug)
!
!
  ENDIF
!
!
!  ----
! END of loop to determine basis parameters
!  -------------------------------------------
!
! in the FIRST segment (nsec=1) only, pev2dmax and potmax
! are set to be exponential (or hyperbolic) functions of rho
!
!     IF(nsec==1)THEN
!     pev2dmax=ap*exp(-alfap*rho)
!     pev2dmax=ap/(rho**alfap)
!     potmax=fmpotma*pev2dmax
!     ENDIF
!  ----------
!
! in the FIRST segment (nsec=1) ONLY, pev2dmax and potmax are
! set to depend on RHO. The first segment is broken into nrhos1-1
! sectors. In each, pev2dmax and potmax are LINEAR functions of RHO.
!
!     IF(nsec/=1) GOTO 90
!
IF(ii==1)THEN
   nrhos1=nrhos1u
   rhos1(nrhos1)=rhos1u(nrhos1u)
   DO nseg=1,nrhos1u-1
      aseg(nseg)=asegu(nseg)
      bseg(nseg)=bsegu(nseg)
      rhos1(nseg)=rhos1u(nseg)
      pevms1(nseg)=pevms1u(nseg)
   ENDDO
ENDIF
!
IF( rho<=rhos1(nrhos1) )THEN
   DO nrpt=2,nrhos1
      IF(rhos1(nrpt)>=rho)THEN
         nseg=nrpt-1
         pev2dmax=aseg(nseg)*rho + bseg(nseg)
         potmax=fmpotma*pev2dmax
         EXIT
       ENDIF
   ENDDO
ENDIF
!
WRITE(Out_Unit,'(/2x,"!! nsec=",i2,2x,"nrho=",i4,2x,"rho=",f15.10," !!")') nsec,nrho,rho
WRITE(Out_Unit,*) 'ntheta,nchi,rhoend = ',ntheta,nchi,rhoend
!
!     IF(nsec==1)THEN
IF( rho>rhos1(nrhos1) )THEN
   WRITE(Out_Unit,*) 'nseg,aseg,bseg = ',nseg,aseg(nseg),bseg(nseg)
ENDIF
!
WRITE(Out_Unit,*) 'potmax,pev2dmax = ',potmax,pev2dmax
WRITE(Out_Unit,'(/)')
!  -----------------------------------------------------
!
WRITE(Out_Unit,'(/2x,"@@ isec=",i2,2x,"nsec=",i2," @@")') isec,nsec
WRITE(Out_Unit,'(2x,"nrho=",i4,2x," rho=",f15.10/)') nrho,rho
!  -------------
nevtot=0
once=.true.
!
!     loop over values of chi and for each angle, calculate a 1-d
!     basis in theta
!
DO itheta=1,ntheta
   theta=acos(pttheta(itheta))
   CALL main1d(rho,potmax,nchi,ltheta,nch,theta,peigv,itheta, debug,h3sys,ii)
   !
   !     at each angle, determine how many eigenvectors are to be
   !     included in the full basis by discarding those
   !     with energy greater than pev2dmax
   !
   DO i=1,nc2tot(itheta)
      IF(e1d(i,itheta)<=pev2dmax)nev(itheta)=nev(itheta)+1
   ENDDO
   nevtot=nevtot+nev(itheta)
   IF(nevtot>nt2max .AND. once)THEN
      WRITE(Out_Unit,*)' ithetmax = ',itheta
      once = .false.
   ENDIF
ENDDO
!
IF(nevtot>nt2max)THEN
   WRITE(Msg_Unit,*)' nevtot2d = ',nevtot,' greater than nt2max', nt2max
   WRITE(Out_Unit,*)' nevtot2d = ',nevtot,' greater than nt2max', nt2max
   STOP "DVR"
ENDIF
!
WRITE(Out_Unit,'(2x,"nrho=",i4,2x,"rho=",f15.10,2x,"nevtot=",i4)') nrho,rho,nevtot
!
!     chctrans calculates the full hamiltonian by adding
!     in the coupling between theta angles. Then it diagonalizes
!     the full hamiltonian htot giving eigenvalues
!     in etot and eigenvectors in cctot
!
CALL chctrans(rho,ntheta,nchi,peigv,debug,ltheta,ngood2)
!
!  ---------------------------------------------------------------------
!  calculate the average of ngood1 and ngood2 surface function
!   eigenvalues
!  ---------------------------------------------------------------------
!
eigenav = 0.0d0
DO i=1,ngood1
  eigenav = eigenav + e2d(i)
ENDDO
eigenav1 = eigenav/ngood1
DO i=ngood1+1,ngood2
  eigenav = eigenav + e2d(i)
ENDDO
eigenav2 = eigenav/ngood2
WRITE(Out_Unit,'("eigenav = ",e16.8," (eV), ngood = ",i4)') eigenav1,ngood1
WRITE(Out_Unit,'("eigenav = ",e16.8," (eV), ngood = ",i4)') eigenav2,ngood2
!
tdum = 0.d0
deltat = tarray(1) - tthen
WRITE(Out_Unit,*) ' '
WRITE(Out_Unit,*) 'time for surface functions = ',deltat,' sec'
!
IF(psurfonl) GOTO 510
!
krho=3
!
IF(ii==1)THEN
   !  upper surface
   !  This set is for unformatted WRITEs to Pmatrx.
   !WRITE(Out_Unit,*)'rho,nsfunc,krho,rhos,mid_zero',rho,nsfunc,krho,rhos,mid_zero
   WRITE(dvr7) rho,nsfunc,krho,rhos,mid_zero
   WRITE(dvr7) (e2d(i)/autoev, i=1,nsfunc)
   !  ---------------------------------------------------------------------
   !  WRITE eigenvalues to sfelevl(u) in FORMAT for plotting.
   !  ---------------------------------------------------------------------
   WRITE(dvr23,'(/," $readx",5x,"distance in bohr")')
   WRITE(dvr23,'(1x,i5)')ione
   WRITE(dvr23,'(1x,f10.5)')rho
   WRITE(dvr23,'(" $ready",5x,"energy in ev.")')
   WRITE(dvr23,'(1x,i5)')nsfunc
   WRITE(dvr23,'(1x,e14.7)')(e2d(i),i=1,nsfunc)
   WRITE(dvr23,'(" $end")')
   OPEN(Unit=DVR_Elevels_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'GraphicsOut/Sfelevl_DVR.csv',Form='formatted',status='unknown')
   WRITE(DVR_Elevels_Unit,'(1x,e14.7,10000(",",e14.7))')rho,(e2d(ibasis),ibasis=1,nsfunc)
ENDIF
!
WRITE(dvr14,*) rho,nevtot

IF(matelem)THEN
   GOTO 1010
ELSE
   GOTO 510
ENDIF
!
 1010 CONTINUE
!
tdum = 0.d0
tthen = tarray(1)
!
!  -------------------------------------------------------
! beginning of the section which ultimately calculates the
! potential coupling matrix in the surface function basis,
! for three rho values within the sector centered at rhoz
!  --------------------------------------------------------
!
DO irho=1,3
IF(irho/=2)THEN
WRITE(Out_Unit,'(/2x,"rhos(",i1,")/au=",f17.13,2x,"rhoz/au=",f17.13/)') irho,rhos(irho),rhoz
!
! loop below calculates potential coupling matrix (eq.168 of pack and
! parker paper) in the primitive 2d dvr basis. elements are stored
! in the array udvr.
!
DO itheta=1,ntheta
   theta=acos(pttheta(itheta))
   nctem=nc2tot(itheta)
   DO i=1,nc2tot(itheta)
      ichi=nptchi(i,itheta)
      chii=acos(ptchi(ichi))
      CALL pptorr(rhoz,theta/Two,chii,rz,dz,chi1z)
      CALL rds2r3(rz,dz,chi1z,dist,nch)
      IF(ii==1)THEN
         PES_Name=pesnamu
         CALL surface(vrhoz,dist)
      ENDIF
      CALL pptorr(rhos(irho),theta/Two,chii,r,d,chi1)
      CALL rds2r3(r,d,chi1,dist,nch)
      IF(ii==1)THEN
         PES_Name=pesnamu
         CALL surface(vrhos,dist)
      ENDIF
      udvr(itheta,i)=autoev*(vrhos-(rhoz/rhos(irho))**2*vrhoz)
   ENDDO
ENDDO
!
! udvr is now transformed (block by block) to the ray eigenvector
! basis, and stored in the urev matrix.
!
DO i=1,nt2max
   DO j=1,nt2max
      urev(i,j)=0.0d0
   ENDDO
ENDDO
!
nevsum=0
!
DO itheta=1,ntheta
   IF(nev(itheta)/=0)THEN
      !
      ! loop over theta blocks
      !
      DO i=1,nc2tot(itheta)
         DO j=1,nev(itheta)
            ucmat(i,j)=udvr(itheta,i)*cc1dtot(i,j,itheta)
         ENDDO
      ENDDO
      DO i=1,nev(itheta)
         DO j=1,nc2tot(itheta)
            ctmat(i,j)=cc1dtot(j,i,itheta)
         ENDDO
      ENDDO
      m=nev(itheta)
      n=nc2tot(itheta)
      l=nev(itheta)
      CALL dgemm('N','N',m ,n ,l ,1.d0 ,ctmat  ,nchimax  ,ucmat  ,nchimax  ,0.d0 ,ctuc  ,nchimax)
      !
      ! ctuc is the theta-th block of the rev potential matrix urev.
      !
      DO i=1,nev(itheta)
         DO j=1,nev(itheta)
            urev(nevsum+i,nevsum+j)=ctuc(i,j)
         ENDDO
      ENDDO
   ENDIF
   nevsum=nevsum+nev(itheta)
ENDDO
!
! END of loop for the formation of the urev matrix
!
! final potential coupling matrix in the surface function basis
! is calculated now. two matrix multiplications are required.
!
m=nevtot
n=nevtot
l=nsfunc
PRINT*,'nevtot = ', nevtot, 'nsfunc = ',nsfunc
CALL dgemm('N','N',m ,n ,l ,1.d0 ,urev   ,nt2max   ,h2d    ,nt2max   ,0.d0 ,urec2 ,nt2max )
!
DO i=1,nevtot
   DO j=1,nevtot
      urev(i,j)=0.0d0
   ENDDO
ENDDO
!
m=nsfunc
n=nevtot
l=nsfunc
CALL dgemm('T','N',m ,l ,n ,1.d0 ,h2d    ,nt2max   ,urec2  ,nt2max   ,0.d0 ,urev  ,nt2max )
!
DO i=1,nsfunc
   DO j=1,nsfunc
      urev(i,j)=urev(i,j)/autoev
   ENDDO
ENDDO
ENDIF
!  ---
! final potential coupling matrix in the surface function basis
! has been calculated and stored in the urev array.
!  ---
IF(irho==2)THEN
   DO i=1,nsfunc
      DO j=1,nsfunc
         urev(i,j)=0.0d0
      ENDDO
   ENDDO
ENDIF
!
WRITE(Out_Unit,'(2x,"potential matrix elements"/)')
!
!  WRITE a few of the matrix elements on 6 IF desired.
IF(medium)THEN
   DO i=1,min(nsfunc,7)
      WRITE(Out_Unit,'(1x,7ES14.6)') (urev(i,j), j=1,min(nsfunc,7))
   ENDDO
   WRITE(Out_Unit,'(2x,"potential matrix elements"/)')
ENDIF
!
IF(irho/=2)THEN
   DO j=1,nsfunc
      IF(ii==1)THEN
         !
         !  This is for an unformatted WRITE.
         !
         WRITE(dvr7) (urev(i,j), i=j,nsfunc)
      ENDIF
   ENDDO
ENDIF
ENDDO
!
tdum = 0.d0
deltat = tarray(1) - tthen
WRITE(Out_Unit,*) 'time for matrix elements = ',deltat,' sec'
!
!  Flush the buffer for units 7, 14, and 23.
!  (8,24)
!
IF(ii==1)THEN
  CALL flush(dvr7)
  CALL flush(dvr14)
  CALL flush(dvr23)
ELSE
  CALL Flush(dvr8)
  CALL FLUSH(dvr24)
ENDIF
!
!  --------------
! END of section calculating potential coupling matrix
! in the surface function basis
!  ---------------------------
!
 510  CONTINUE
!
!
IF(ii==1)THEN
   REWIND nunita
   WRITE(nunita)nrho,rho,nsfunc,ntheta,nchi
!  WRITE(nunita,*)nrho,rho,nsfunc,ntheta,nchi
ENDIF
!
!
!  ------------------------------------------------------
! surface functions (nsfunc of them) are now, one by one,
! transformed to primitive dvr basis
!  ------------------------------------------------------
!
IF(nrho /= 1 .AND. (.NOT. psurfonl) )THEN
  tdum = 0.d0
  tthen = tarray(1)
ENDIF
!
nthnch=ntheta*nchi
!
DO ksf=1,nsfunc
   DO i=1,ndvptmax
      csurf1(i,1)=0.0d0
   ENDDO
   nreth=0
   nchith=0
   DO itheta=1,ntheta
      IF(nev(itheta)==0)THEN
         DO m=1,nchi
            csurf1(nchith+m,1)=0.0d0
         ENDDO
         GOTO 670
      ENDIF
      DO m=1,nchi
         sumc=0.0d0
         DO ip=1,nev(itheta)
            sumc=sumc+cc1d(m,ip,itheta)*h2d(nreth+ip,ksf)
         ENDDO
         csurf1(nchith+m,1)=sumc
      ENDDO
!
 670  CONTINUE
      nreth=nreth+nev(itheta)
      nchith=nchith+nchi
   ENDDO
   !  ---
   ! ksf-th surface eigenvector in the primitive dvr basis is
   ! written out to unit 12 or unit 13, depending on the value
   ! of nunita, which changes from one rho to another
   !  ---
   !
   IF(ii==1)THEN
      WRITE(nunita) (csurf1(i,1), i=1,nthnch)
      ! WRITE(nunita,*) (csurf1(i,1), i=1,nthnch)
      !IF(nrho==1.AND.ksf==1)THEN
         !ALLOCATE(ThetaVal(nthnch),Chivals(nthnch))
         !index=0
         !DO itheta=1,ntheta
            !DO m=1,nchi
              ! index=index+1
              ! ThetaVal(index)=ACOS((pttheta(itheta)+1.d0)/2.d0)
              ! ChiVals(index)=ACOS(ptchi(m)/2)
            !ENDDO
         !ENDDO
         !SFTYPE="Temp01"
         !CALL SFPlots(nthnch, ThetaVal, ChiVals, CSurf1, SFType)
         !Write(msg_unit,*)'ThetaVal=',ThetaVal
         !Write(msg_unit,*)'ChiVals=',ChiVals
         !DEALLOCATE(ThetaVal,ChiVals)
         !STOP "TEMPP"
      !ENDIF
   ENDIF
ENDDO
!
!  --------------------------------------------------------------------
! WRITE out the transformtion matrix for theta and for chi
!  --------------------------------------------------------------------
REWIND nunitta
WRITE(nunitta) rho
WRITE(nunitta) ((ttheta(i,j),i=1,ntheta),j=1,ntheta)
WRITE(nunitta) ((tchi(i,j),i=1,nchi),j=1,nchi)
!  ---------------------------------------------------------------------
! END of section transforming surface eigenfunctions
! to the primitive dvr basis
!  ---------------------------------------------------------------------
!
IF(psurfonl) GOTO 1940
IF(.NOT.matelem) GOTO 1940
!  --------------------------------------------------------------------
! overlap matrix not calculated at first rho
!  --------------------------------------------------------------------
IF(nrho==1) GOTO 1940
!  --------------------------------------------------------------------
!  calculate Ovrlp matrix between surface functions at the current and
!  previous rho.
!  --------------------------------------------------------------------
IF(nchi == nchip .AND. ntheta == nthetap.or.nrho==irhost)THEN
!  --------------------------------------------------------------------
!   case of identical dvr orders for the two values of rho.
!  --------------------------------------------------------------------
  IF(ii==1)THEN
    CALL ovrcald (nunitb, nunita, nsfunc, nthnch,nrho, rholast, rho,ii, medium)
  ENDIF
  tdum = 0.d0
  deltat = tarray(1) - tthen
  WRITE(Out_Unit,*) 'time for overlaps (1st)  = ',deltat,' sec'
ELSE
!
!  ---------------------------------------------------------------
! Calculate overlaps matrices for the case where the DVR orders
!  are NOT identical for the two values of rho
! note that this is currently set to be the same for upper
! and lower surfaces.  Not complete.
!  ---------------------------------------------------------------
  IF(ii==1)THEN
     CALL ovrcal2 (nunitb, nunita, nsfunc, ntheta, nthetap,nchi, nchip, nrho, rholast, rho)
  ENDIF
!
!
!
  tdum = 0.d0
  deltat = tarray(1) - tthen
  WRITE(Out_Unit,*) 'time for overlaps (2nd)  = ',deltat,' sec'
!
ENDIF
!
!  at last rho value, WRITE one more line to Ovrlp for the propagator.
!
IF(nrho==nrhomx)THEN
   WRITE(dvr18) nrho, rho, rhonext                           !Be Careful
ENDIF
!
!  Flush the buffer for unit 19 or 20 (Ovrlp).
!
IF(ii==1)THEN
  CALL flush(dvr18)
ENDIF
!
!  --------------------------------------------------------------------
! END of overlap matrix formation
!  --------------------------------------------------------------------
!
 1940 CONTINUE
!
!  --------------------------------------------------------------------
!  decide whether to DO dvr to fbr transformation.
!  the fbr is needed to transform to the asymptotic or abm functions.
!  --------------------------------------------------------------------
!
!IF(rho > rhofbr1 .AND.rho < rhofbr2.AND.ii==1)THEN
IF(nrho==nrhomx)THEN
   !
   !  --------------------------------------------------------------------
   !  begin the DVR to FBR transformation
   !  --------------------------------------------------------------------
   !
   tdum = 0.d0
   tthen = tarray(1)
   !
   !  Both theta and chi DVR transformation matrices, ttheta and tchi
   !  are transposed, to conform to the usual notation (first index
   !  labels the basis function, second labels the angle).
   !
   WRITE(Out_Unit,'(/2x,"DVR to FBR transformation executing", 4x,"rho = ",f15.10/)') rho
   DO i=1,ntheta
      DO j=1,ntheta
         tthetat(i,j) = ttheta(j,i)
      ENDDO
   ENDDO
   DO i=1,nchi
      DO j=1,nchi
         tchit(i,j) = tchi(j,i)
      ENDDO
   ENDDO
   !
   !  This is for an unformatted WRITE,
   !
   ! redone by m.b. to try to avoid problems, 9/91
   !
   WRITE(dvr22) nrho,rho,ntheta,nchi
   WRITE(dvr22) (nev(jj),jj=1,ntheta)
   DO jj=1,ntheta
     WRITE(dvr22) (tthetat(jj,kk),kk=1,ntheta)
   ENDDO
   !
   !  Coefficients of the surface eigenvectors (nsfunc of them) in
   !  the chi DVR basis are, one by one, transformed to the coefficients
   !  in the chi FBR basis
   !
   REWIND nunita
   READ(nunita) nrhdum,rhodum,nsfdum,nthdum,nchdum
   !
   DO ksf=1,nsfunc
      READ(nunita) (csurf1(k,1),k=1,nthnch)
      nchith = 0
      DO itheta=1,ntheta
         IF(nev(itheta) == 0)THEN
            DO m=1,nchi
               csurfbr(nchith + m) = 0.0d0
            ENDDO
         ELSE
            DO m=1,nchi
               sumc = 0.0d0
               DO k=1,nchi
                  sumc = sumc + tchit(m,k)*csurf1(nchith + k,1)
               ENDDO
               csurfbr(nchith + m) = sumc
            ENDDO
         ENDIF
         nchith = nchith + nchi
      ENDDO
      WRITE(dvr22) (csurfbr(i),i=1,nthnch)
   ENDDO
      
   !
   CALL flush(dvr22)
   !
   tdum = 0.d0
   deltat = tarray(1) - tthen
   WRITE(Out_Unit,*) 'time for FBR transform   = ',deltat,' sec'
   !
   CLOSE(dvr22)
   !
   !  -------------------------------------------------------------------
   !  this is the ENDIF that ends the dvr to fbr transformation section
   !  -------------------------------------------------------------------
   
   !IF(nrho==20)THEN
   Ntheta_plt=31
   NChi_plt=31
   ALLOCATE(ThetaVal(Ntheta_plt*NChi_plt))
   ALLOCATE(Chivals(Ntheta_plt*NChi_plt))
   ALLOCATE(Phi_Plt(Ntheta_plt*NChi_plt))
   index=0
   DO itheta=1,Ntheta_plt
      DO m=1,NChi_plt
         index=index+1
         ThetaVal(index)=(itheta-1)*Pi/(2.d0*(Ntheta_plt-1))
         ChiVals(index)=(m-1)*2.d0*Pi/(NChi_plt-1)
      ENDDO
   ENDDO
   SFTYPE="Temp01"
   CALL Phi_DVR(ChiVals, NChi, ThetaVal, NTheta, Phi_Plt, Ntheta_plt*NChi_plt, 2, rho)
   CALL SFPlots(Ntheta_plt*NChi_plt, ThetaVal, ChiVals, Phi_Plt, SFType)
   DEALLOCATE(ThetaVal,ChiVals,Phi_Plt)
   !STOP "TEMPP"
!ENDIF


   isfn = isfn + 1
   !  This set is for unformatted WRITEs.
   OPEN(Unit=dvr22,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/SfunFBR.bin',status='unknown', form='unformatted')

ENDIF
!  -------------------------------------------------------------------
isec=0
psurfonl=.false.
!
IF(rho==rhofinal) GOTO 599
!
! Store the transformation matrices for theta and chi into new
!  arrays for use in computing overlaps for non-identical DVR orders.
!
nthetap = ntheta
nchip = nchi
!
DO i=1,ntheta
   DO j=1,ntheta
      tthetap(i,j) = ttheta(i,j)
   ENDDO
ENDDO
DO i=1,nchi
   DO j=1,nchi
      tchip(i,j) = tchi(i,j)
   ENDDO
ENDDO
!
! -------------------------------------------------------------------
! END loop over upper and lower PE surfaces (ii=1, nsurf)
!  ------------------------------------------------------------------

 1960 CONTINUE

!
!  -------------------------------------------------------------------
!  END of dipole matrix element calculation
!  -------------------------------------------------------------------
!
!  *******************************************************************
! END OF LOOP OVER RHO VALUES
!  *******************************************************************
!
ENDDO
!
!  -------------------------------------------------------------------
 599  CONTINUE
!
CLOSE(dvr7)
CLOSE(dvr18)
CLOSE(dvr23)
CLOSE(Out_Unit)
OPEN(Unit=Out_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'Output/BasisForScatt.txt',Form='FORMATTED',ACCESS='APPEND')
WRITE(Out_Unit,*)'Completed DVR'
WRITE(Out_Unit,*)

RETURN
ENDSUBROUTINE DVR
