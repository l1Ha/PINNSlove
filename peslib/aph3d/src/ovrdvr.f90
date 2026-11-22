SUBROUTINE ovrdvr (Rho_DVR, Rho_ABM, Nang_ABM, nbasis, xchanl, nangch, Phi_DVR, Phi_ABM, &
                   vect, Overlap_DVR2ABM, amatrx, nbasiss,                               &
                   neigmin, jrot,  Theta_ABM, Chi_ABM, Weight_ABM)
USE FileUnits_Module
USE Parms_Module
USE TotJ_Module
USE FBR_Module
USE Testing_Module
USE InputFile_Module
!

!            P U R P O S E   O F    S U B R O U T I N E
! Determines dvr to abm overlap matrix elements
! using abm quadratures and angles.


!            I N P U T    A R G U M E N T S
!   Rho_DVR   previous hyperradius.
!   Rho_ABM   current hyperradius
!   Nang_ABM  number of ABM quadrature angles.
!   nbasis    number of raw primitive ABM basis functions.
!   nbasiss   number of symmetry-adapted primitive ABM basis functions.
!   narran    number of arrangement channels.
!   xchanl    an array giving the channel number of each basis function.
!   nangch    an array containing the starting position of the
!             angles for each of the arrangement channels.
!             nangch(narran+1) is equal to Nang_ABM+1.
!   Phi_DVR   primitive basis at Rho_DVR at ABM quadrature pts.
!   Phi_ABM   primitive basis at Rho_ABM at ABM quadrature pts.
!   vect      coeffs of present surface functions.

!            O U T P U T    A R G U M E N T S
!   Overlap_DVR2ABM        overlap matrix elements.


IMPLICIT NONE
CHARACTER(LEN=200) filename, prefix

!            I N T E G E R S
INTEGER Nang_ABM, ibasis, jbasis, nbasis, iang, iangst1, iangend1, jpar, neigmin, lnblnk
INTEGER nbasiss, nmode,  nprimabm, nprimdvr, nabm, ndvr, isav, IErr
INTEGER xchanl(nbasiss), nangch(narran+1), jrot(nbasiss)

!            R E A L S
REAL(Kind=WP_Kind) overlap, Rho_DVR, Rho_ABM, zero, half, two, xnorm, xfac
REAL(Kind=WP_Kind) Phi_DVR(Nang_ABM, nbasiss), Overlap_DVR2ABM(nbasiss, nbasiss), Phi_ABM(Nang_ABM, nbasiss)
REAL(Kind=WP_Kind) vect(nbasiss, nbasiss), amatrx(nbasiss,nbasiss)
REAL(Kind=WP_Kind) Theta_ABM(Nang_ABM), Chi_ABM(Nang_ABM), Weight_ABM(Nang_ABM)
REAL(Kind=WP_Kind) diagonal(1000)

!            E X T E R N A L S
EXTERNAL popt, mxoutf, sgemm_junk, rdwrsurb

!            P A R A M E T E R S
PARAMETER (zero=0.0d0, half=0.5d0, two=2.0d0)
!-----------------------------------------------------------------------
! Determine printing options.
!-----------------------------------------------------------------------
LOGICAL little, medium, full, there
INTEGER ithcall, ithsub
DATA ithcall /0/, ithsub /0/, isav /-1/
DATA little /.false./, medium /.false./, full /.false./
CALL popt ('ovrdvr  ', little, medium, full, ithcall, ithsub)

xnorm = sqrt(half)
xfac = sqrt(two)

INQUIRE(file=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile, exist = there)
IF(there)THEN
   OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile,status='old')
   READ(In_Unit,NML=testing, IOSTAT=IERR)
   IF(IERR==-1)THEN
      test='nul'
   ENDIF
   CLOSE(Unit=In_Unit)
   WRITE(Out_Unit,NML=testing)
ELSE
   STOP 'Stopping in DVR_to_ABM'
ENDIF

IF(test=='ABM')THEN
   WRITE(Msg_Unit,*)'Testing ABM Surface Functions'
   WRITE(Out_Unit,*)'Testing ABM Surface Functions'
   WRITE(Msg_Unit,*)'Warning!  Not normalized for this case!'
   WRITE(Out_Unit,*)'Warning!  Not normalized for this case!'
ELSEIF(test=='DVR')THEN
   WRITE(Msg_Unit,*)'Testing DVR Surface Functions'
   WRITE(Out_Unit,*)'Testing DVR Surface Functions'
   WRITE(Msg_Unit,*)'Warning!  Not normalized for this case!'
   WRITE(Out_Unit,*)'Warning!  Not normalized for this case!'
ELSE
   WRITE(Msg_Unit,*)'Calculating DVR_TO_ABM transformation'
   WRITE(Out_Unit,*)'Calculating DVR_TO_ABM transformation'
ENDIF
!---------------------------------------------------------------
!  READ the number of fbr(dvr) functions to be used.
!---------------------------------------------------------------

INQUIRE(file=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile,exist=there)
IF(there)THEN
  OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile,status='old')
  READ(In_Unit,NML=fbr, IOSTAT=IERR)
   IF(IERR==-1)THEN
      CALL NameList_Default(Out_Unit, 'FBR')
      NSfunc=105 !TmpModGregParker2
      lmax=20
      mmax=40
      nsym=1
   ELSEIF(IERR/=0)THEN
      WRITE(Out_Unit,*)'ERROR with Namelist FBR'
      WRITE(Msg_Unit,*)'ERROR with Namelist FBR'
      STOP 'OvrDVR: FBR'
   ENDIF
  nmode = nsfunc
  CLOSE(Unit=In_Unit)
  WRITE(Out_Unit,NML=fbr)
ELSE
  WRITE(Msg_Unit,*)'Error InputFile does not exist (ovrdvr):',InputDIR(1:LEN(TRIM(InputDIR)))//InputFile
  STOP 'ovrdvr'
ENDIF

WRITE(Out_Unit,*)'nangch=',nangch,' Nang_ABM=',Nang_ABM

!----------------------------------------------------------------
!  construct the fbr(dvr) functions at the abm quadrature points.
!----------------------------------------------------------------

CALL intdvr (Nang_ABM, Chi_ABM, Theta_ABM, Phi_DVR, nmode, Rho_DVR)
CALL weight_unset (Nang_ABM, nprimabm, Phi_ABM, Weight_ABM)

!-------------------------------------------------------------------
! The following only used to isolate inaccuracies in ABM or DVR basis.
!-------------------------------------------------------------------
IF(test=='ABM')THEN
   WRITE(Msg_Unit,*)'Both Phi_DVR and Phi_ABM are ABM s.f.'
   WRITE(Out_Unit,*)'Both Phi_DVR and Phi_ABM are ABM s.f.'
   CALL vcopy(Nang_ABM*NBasiss, Phi_DVR, Phi_ABM)
   nabm = neigmin
   ndvr = neigmin
   nprimdvr = nbasiss
   nprimabm = nbasiss
ELSEIF(test=='DVR')THEN
   WRITE(Msg_Unit,*)'Both Phi_DVR and Phi_ABM are DVR s.f.'
   WRITE(Out_Unit,*)'Both Phi_DVR and Phi_ABM are DVR s.f.'
   CALL vcopy(Nang_ABM*NBasiss, Phi_ABM, Phi_DVR)
   nabm = nmode
   ndvr = nmode
   nprimdvr = nmode
   nprimabm = nmode
ELSE
   WRITE(Msg_Unit,*)'dvr to abm case'
   WRITE(Out_Unit,*)'dvr to abm case'
   nabm = neigmin
   ndvr = nmode
   nprimdvr = nmode
   nprimabm = nbasiss
ENDIF

IF(ndvr>NBasiss.or.nabm>NBasiss.or.NEigmin>NBasiss.or.NMode>NBasiss)THEN  
   WRITE(Out_Unit,*)'ndvr>NBasiss.or.nabm>NBasiss.or.NEigmin>NBasiss.or.NMode>NBasiss'
   WRITE(Out_Unit,*)ndvr,nabm,NEigmin,NMode,NBasiss
   WRITE(Msg_Unit,*)'ndvr>NBasiss.or.nabm>NBasiss.or.NEigmin>NBasiss.or.NMode>NBasiss'
   WRITE(Msg_Unit,*)ndvr,nabm,NEigmin,NMode,NBasiss
   STOP "Stopping in OvrDVR"
ENDIF

WRITE(Out_Unit,*)'Number of ABM Surface Functions =',nabm
WRITE(Out_Unit,*)'Number of DVR Surface Functions =',ndvr
WRITE(Out_Unit,*)'Number of Primitive ABM Basis Functions=',nprimabm
WRITE(Out_Unit,*)'Number of Primitive DVR Basis Functions=',nprimdvr

CALL SFPlots(Nang_ABM, Theta_ABM, Chi_ABM, Phi_DVR, "DVR_AA")
CALL SFPlots(Nang_ABM, Theta_ABM, Chi_ABM, Phi_ABM, "ABM_AA")

!
! Multiply the ABM and DVR surface functions by the weight.
!
CALL weight_set (Nang_ABM, nprimdvr, Phi_DVR, Weight_ABM)
CALL weight_set (Nang_ABM, nprimabm, Phi_ABM, Weight_ABM)

! This will print weighted surface functions.
medium=.true.
IF(medium)THEN
   WRITE(Out_Unit,*)'chanl, iang, Phi_DVR, Phi_ABM'
   iangst1 = nangch(xchanl(1))
   iangend1 = nangch(xchanl(1)+1)-1
   WRITE(Out_Unit,*)'iangst1,iangend1',iangst1,iangend1
   WRITE(Out_Unit,'(2x,A,3x,A,4x,A,5x,A)')"xchanl", "iang", "Phi_DVR", "Phi_ABM"
   DO iang = iangst1, iangend1
      IF(ABS(Phi_DVR(iang,1))>0.05.or. ABS(Phi_ABM(iang,1))>0.05)THEN
         WRITE(Out_Unit,'(2I5,2ES15.7)')xchanl(1), iang, Phi_DVR(iang,1),Phi_ABM(iang,1)
      ENDIF
   ENDDO
ENDIF


!
!  --------------------------------------------------------------
!  evaluate overlaps of dvr surface functions and primitive abm basis.
!  --------------------------------------------------------------

WRITE(Msg_Unit,*)'Calculate overlap integrals'
WRITE(Out_Unit,*)'Calculate overlap integrals'
WRITE(Out_Unit,*)'Symmetry=',symmetry
IF(symmetry)THEN
   IF(jeven)THEN
      jpar = 0
   ELSE
      jpar = 1
   ENDIF
ENDIF
DO ibasis = 1, nprimabm             !  start loop over abm primitive basis
   iangst1 =  nangch(xchanl(ibasis))
   iangend1 = nangch(xchanl(ibasis)+1)-1
   IF(iangst1/=isav.or.ibasis==nprimabm)WRITE(Out_Unit,*)'ibasis=',ibasis,iangst1,iangend1,xchanl(ibasis)
   isav=iangst1
   DO jbasis =1,nprimdvr            !  start loop over dvr surface functions
      overlap = 0.0d0
      DO iang = iangst1, iangend1   !  quadrature loop
         overlap = overlap + Phi_DVR(iang, jbasis)* Phi_ABM(iang,ibasis)
      ENDDO
!
!  The following makes use of the fact that the dvr surface functions
!  and the abm primitives in channel 1 have the full symmetry, and the
!  quadrature in channel 1 automatically covers half the space and doubles
!  the result, so no renormalization is needed.
!  the quadratures done in channel 2 cover the whole space but only use
!  one of the pair of primitives necessary to have proper symmetry.
!  They must be multiplied by 2/sqrt2.  R. T Pack, Feb. 1993.
!
      IF(symmetry.AND.xchanl(ibasis)==2) overlap=xfac*overlap
      amatrx(jbasis, ibasis) = overlap
   ENDDO
ENDDO
WRITE(Msg_Unit,*)'Overlap integrals have been calculated'
WRITE(Out_Unit,*)'Overlap integrals have been calculated'

WRITE(Out_Unit,*)'Overlaps of Primitive Basis in ovrdvr'
WRITE(Out_Unit,*)'Matrix dimensions:',nprimdvr,'  by',nprimabm
CALL mxoutf (amatrx, nprimdvr, nprimabm, NBasiss, NBasiss)

!--------------------------------------------------------------------
! Transform primitive overlaps integrals to actual overlaps.
!--------------------------------------------------------------------

IF(test/='DVR')THEN
   WRITE(Out_Unit,*)'Coefficients of the primitive basis'
   CALL mxoutf(vect, nprimdvr, nabm, NBasiss, Nabm)
ENDIF

IF(test/='DVR')THEN
   WRITE(Msg_Unit,*)'Transform primitive overlaps to actual overlaps'
   WRITE(Out_Unit,*)'Transform primitive overlaps to actual overlaps'
   CALL DGEMM('N','N',nprimdvr,nabm,nprimabm,1.d0,amatrx,nprimabm,vect,NBasiss,0.d0,Overlap_DVR2ABM,NBasiss)
ELSEIF(test=='DVR')THEN
   Overlap_DVR2ABM=0.d0
   CALL vcopy(Nbasiss*NBasiss, Overlap_DVR2ABM, amatrx)
ENDIF

IF(test=='ABM')THEN
   CALL DGEMM('T','N',nabm,nabm,nprimabm,1.d0,vect,nabm,Overlap_DVR2ABM,nprimabm,0.d0,amatrx,NBasiss)
   Overlap_DVR2ABM=0.d0
   CALL vcopy(NBasiss*NBasiss, Overlap_DVR2ABM, amatrx)
ENDIF

!-----------------------------------------------------------------------
! IF desired WRITE out the overlap matrix.
!-----------------------------------------------------------------------
WRITE(Out_Unit,*)'routine ovrdvr'
WRITE(Out_Unit,*)'ndvr = ', ndvr
WRITE(Out_Unit,*)'nabm = ', nabm

WRITE(Out_Unit,*)'Overlaps Actual basis'

WRITE(Out_Unit,*) 'Part of the overlap matrix at Rho_DVR, Rho_ABM = ', Rho_DVR, Rho_ABM
WRITE(Out_Unit,*)'Matrix dimensions:',ndvr,'  by',nabm
CALL mxoutf(Overlap_DVR2ABM, ndvr, nabm, NBasiss, NBasiss)


!  --------------------------------------------------------------
!  WRITE overlap matrix on unit Ovr_Unit 'Ovrlp'
!  --------------------------------------------------------------
WRITE(Ovr_Unit) 0, Rho_DVR, Rho_ABM
WRITE(Ovr_Unit) nabm, Rho_DVR, Rho_ABM
DO jbasis=1,nabm
   WRITE(Ovr_Unit)(Overlap_DVR2ABM(ibasis,jbasis),ibasis=1,ndvr)
ENDDO
WRITE(Out_Unit,'(A)')'Diagonal Elements of Overlap basis'
DO ibasis = 1, min(ndvr,nabm)
   diagonal(ibasis)=Overlap_DVR2ABM(ibasis,ibasis)
ENDDO
WRITE(Out_Unit,'(10ES11.4)')(diagonal(ibasis),ibasis=1,min(ndvr,nabm))

!--------------------------------------------------------------------
!  construct Overlap_DVR2ABM*Overlap_DVR2ABM(trans) to test completeness of new basis.
!--------------------------------------------------------------------
WRITE(Out_Unit,*) 'new basis completeness test.  Rho_DVR, Rho_ABM= ',Rho_DVR, Rho_ABM
CALL vsets(NBasiss*NBasiss, amatrx, 1, 0.d0)
CALL DGEMM('N','T',ndvr,ndvr,nabm,1.d0,Overlap_DVR2ABM,NBasiss,Overlap_DVR2ABM,NBasiss,0.d0,amatrx,NBasiss)
CALL mxoutf (amatrx, ndvr, ndvr, NBasiss, NBasiss)
WRITE(Out_Unit,'(A)')'Diagonals of new basis completeness test'
DO ibasis = 1, ndvr
   diagonal(ibasis)=amatrx(ibasis,ibasis)
ENDDO
WRITE(Out_Unit,'(10ES11.3)')(diagonal(ibasis),ibasis=1,ndvr)
!--------------------------------------------------------------------
!  construct Overlap_DVR2ABM(trans)*Overlap_DVR2ABM to test completeness of old basis.
!--------------------------------------------------------------------
WRITE(Out_Unit,'(A,2ES15.7)') 'old basis completeness test at Rho_DVR, Rho_ABM= ', Rho_DVR, Rho_ABM
CALL vsets(NBasiss*NBasiss, amatrx, 1, 0.d0)
CALL DGEMM('T','N',nabm,ndvr,nabm,1.d0,Overlap_DVR2ABM,NBasiss,Overlap_DVR2ABM,NBasiss,0.d0,amatrx,NBasiss)
CALL mxoutf (amatrx, nabm, nabm, NBasiss, NBasiss)
WRITE(Out_Unit,'(A)')'Diagonals of old basis completeness test'
DO ibasis = 1, nabm
   diagonal(ibasis)=amatrx(ibasis,ibasis)
ENDDO
WRITE(Out_Unit,'(10ES11.3)')(diagonal(ibasis),ibasis=1,nabm)

IF(test=='ABM')THEN
   WRITE(Out_Unit,*)'END of ABM test'
   STOP 'ABM test'
ELSEIF(test=='DVR')THEN
   WRITE(Out_Unit,*)'END of DVR test'
   STOP 'DVR test'
ELSE
   WRITE(Out_Unit,'(A)')'Calculated Matrix between DVR and ABM Surface Functions'
ENDIF

RETURN
END
