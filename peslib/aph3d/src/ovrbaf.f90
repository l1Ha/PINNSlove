SUBROUTINE ovrbaf (rholast, rho, nbasis, narran, vect, s, amatrx, neigmin, megamin, mega, nbasiss)
USE FileUnits_Module
USE Numbers_Module
USE Testing_Module
USE InputFile_Module
!            P U R P O S E   O F    S U B R O U T I N E
! This routine determines the overlap matrix elements since basis
! set is a nonorthogonal basis.
!            I N P U T    A R G U M E N T S
!   rholast  previous hyperradius.
!   rho      current hyperradius
!   nbasis   number of basis functions.
!   narran   number of arrangement channels.
!   vect     coeffs of present surface functions.

!            O U T P U T    A R G U M E N T S
!   s        overlap matrix elements.

IMPLICIT NONE
LOGICAL there
CHARACTER(LEN=200) prefix
INTEGER  ibasis, jbasis, nbasis, narran, neigmin, megamin, mega, irhold, ndvr, nabm
INTEGER  ndimen, lnblnk, nprimabm, nprimdvr, nbasiss
REAL(Kind=WP_Kind) rholast, rho, diagonal(1500)
REAL(Kind=WP_Kind) s(nbasiss, nbasiss), vect(nbasiss, nbasiss), amatrx(nbasiss,nbasiss)

EXTERNAL popt, mxoutd, sgemm_junk, rdwrsurb

!-----------------------------------------------------------------------
! Determine printing options.
!-----------------------------------------------------------------------
LOGICAL little, medium, full
INTEGER ithcall, ithsub
COMMON /rhoold/ irhold
DATA ithcall /0/, ithsub /0/
DATA little /.false./, medium /.false./, full /.false./
CALL popt ('ovrbaf  ', little, medium, full, ithcall, ithsub)
little = .true.

WRITE(Out_Unit,*)'ovrbaf called:'
WRITE(Out_Unit,*)'rholast,rho=',rholast,rho
WRITE(Out_Unit,*)'nbasis=',nbasis

ndimen = nbasiss
INQUIRE(file=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile, exist = there)
IF(there)THEN
   OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile,status='old')
   READ(In_Unit,NML=testing)
   CLOSE(Unit=In_Unit)
   WRITE(Out_Unit,NML=testing)
ELSE
   STOP "OvrBaf"
ENDIF

IF(test=='ABM')THEN
   WRITE(Msg_Unit,*)'Testing ABM Surface Functions'
   WRITE(Out_Unit,*)'Testing ABM Surface Functions'
ELSEIF(test=='DVR')THEN
   WRITE(Msg_Unit,*)'Testing DVR Surface Functions'
   WRITE(Out_Unit,*)'Testing DVR Surface Functions'
ELSE
   WRITE(Msg_Unit,*)'Calculating DVR_TO_ABM transformation'
   WRITE(Out_Unit,*)'Calculating DVR_TO_ABM transformation'
ENDIF

CALL vsets(ndimen*ndimen, s, 1, 0.d0)
WRITE(Out_Unit,*) 'ovrbaf: nbasis=',nbasiss
!----------------------------------------------------------------
! READ in DVR surface functions.
! calculate ABM surface functions at nodal points.
! construct overlap matrix.
!----------------------------------------------------------------
CALL basbaf (nbasis, ndvr, neigmin, vect, amatrx, rholast, irhold, test, nbasiss)

!-------------------------------------------------------------------
! The following is only used to isolate inaccuracies in ABM or DVR basis.
!-------------------------------------------------------------------
IF(test=='ABM')THEN
   WRITE(Msg_Unit,*)'Both Phi_DVR and Phi_ABM are ABM s.f.'
   WRITE(Out_Unit,*)'Both Phi_DVR and Phi_ABM are ABM s.f.'
   nabm = neigmin
   ndvr = neigmin
   nprimdvr = nbasiss
   nprimabm = nbasiss
ELSEIF(test=='DVR')THEN
   WRITE(Msg_Unit,*)'Both Phi_DVR and Phi_ABM are DVR s.f.'
   WRITE(Out_Unit,*)'Both Phi_DVR and Phi_ABM are DVR s.f.'
   nabm = ndvr
   ndvr = ndvr
   nprimdvr = ndvr
   nprimabm = ndvr
ELSE
   nabm = neigmin
   ndvr = ndvr
   nprimdvr = ndvr
   nprimabm = nbasiss
ENDIF

!--------------------------------------------------------------------
! Transform primitive overlaps integrals to actual overlaps.
!--------------------------------------------------------------------
WRITE(Out_Unit,*)'Overlaps of Primitive Basis in ovrbaf'
WRITE(Out_Unit,*)'Matrix dimensions:',nprimdvr,'  by',nprimabm
CALL mxoutf (amatrx,nprimdvr,nprimabm,ndimen,ndimen)

IF(test/='DVR')THEN
   WRITE(Out_Unit,*)'Coefficients of the primative basis'
   CALL mxoutf(vect,nprimabm,nabm,ndimen,ndimen)
ENDIF

IF(test/='DVR')THEN
   WRITE(Msg_Unit,*)'Transform primitive overlaps to actual overlaps'
   WRITE(Out_Unit,*)'Transform primitive overlaps to actual overlaps'
   CALL sgemm_junk(nprimdvr, nprimabm, nabm, amatrx, ndimen, vect, ndimen, s, ndimen, 0 ,1)
ELSEIF(test=='DVR')THEN
   CALL vcopy(ndimen*ndimen,s,amatrx)
ENDIF

IF(test=='ABM')THEN
   CALL sgemm_junk(nabm, nprimabm, nabm, vect, ndimen, s, ndimen, amatrx, ndimen,4,1)
   CALL vcopy(ndimen*ndimen,s,amatrx)
ENDIF

WRITE(Out_Unit,*)'Overlaps Actual basis'
WRITE(Out_Unit,*)'Matrix dimensions:',ndvr,'  by',nabm
CALL mxoutf(s,ndvr,nabm,ndimen,ndimen)


!  --------------------------------------------------------------
!  WRITE overlap matrix on unit Ovr_Unit 'Ovrlp'
!  --------------------------------------------------------------
WRITE(Ovr_Unit) irhold, rholast, rho
WRITE(Ovr_Unit) neigmin, rholast, rho
DO jbasis=1,neigmin
   WRITE(Ovr_Unit)(s(ibasis,jbasis),ibasis=1,neigmin)
ENDDO

WRITE(Out_Unit,'(A)')'Diagonal Elements of Overlap basis'
DO ibasis = 1, min(ndvr,nabm)
   diagonal(ibasis)=s(ibasis,ibasis)
ENDDO
WRITE(Out_Unit,'(10ES11.4)')(diagonal(ibasis),ibasis=1,min(ndvr,nabm))

!-----------------------------------------------------------------------
! IF desired WRITE out the overlap matrix.
!-----------------------------------------------------------------------
WRITE(Out_Unit,*)'routine ovrbaf'
WRITE(Out_Unit,*)'ndvr = ', ndvr
WRITE(Out_Unit,*)'nabm =', nabm

WRITE(Out_Unit,*) 'Part of the overlap matrix at rholast, rho = ', rholast, rho
WRITE(Out_Unit,*)'Matrix dimensions:',ndvr,'  by',nabm
CALL mxoutf (s, ndvr, nabm, ndimen, ndimen)
!--------------------------------------------------------------------
!  construct s*s(trans) to test completeness of new basis.
!--------------------------------------------------------------------
WRITE(Out_Unit,*) 'new basis completeness test.  rholast, rho= ',rholast, rho
CALL vsets(ndimen*ndimen, amatrx, 1, 0.d0)
CALL sgemm_junk(ndvr, nabm, ndvr, s, ndimen,s, ndimen, amatrx, ndimen, 2, 1)
CALL mxoutf (amatrx, ndvr, ndvr, ndimen, ndimen)
WRITE(Out_Unit,*)'Diagonal Contributions'
DO ibasis = 1, ndvr
   diagonal(ibasis)=amatrx(ibasis,ibasis)
ENDDO
!WRITE(Out_Unit,*)(diagonal(ibasis),ibasis=1,ndvr)
!--------------------------------------------------------------------
!  construct s(trans)*s to test completeness of old basis.
!--------------------------------------------------------------------
WRITE(Out_Unit,*) 'old basis completeness test.  rholast, rho= ',rholast, rho
CALL vsets(ndimen*ndimen, amatrx, 1, 0.d0)
CALL sgemm_junk(ndvr, nabm, nabm, s, ndimen,s, ndimen, amatrx, ndimen, 4, 1)
CALL mxoutf (amatrx, ndvr, nabm, ndimen, ndimen)
WRITE(Out_Unit,*)'Diagonal Contributions'
DO ibasis = 1, nabm
   diagonal(ibasis)=amatrx(ibasis,ibasis)
ENDDO
!WRITE(Out_Unit,*)(diagonal(ibasis),ibasis=1,nabm)

IF(test=='ABM')THEN
   WRITE(Out_Unit,*)'END of ABM test'
   STOP 'ABM test'
ELSEIF(test=='DVR')THEN
   WRITE(Out_Unit,*)'END of DVR test'
   STOP 'DVR test'
ELSE
   WRITE(Out_Unit,*)'Calculated Matrix between DVR and ABM ','Surface Functions'
ENDIF

RETURN
END
