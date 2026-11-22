SUBROUTINE OvrFEM (rholast, rho, nang, nbasis, xchanl, nangch, Phi_FEM, &
  Phi_ABM, vect, s, amatrx, nbasiss, neigmin, jrot, ThetAPH, chiaph, weight)
USE Numbers_Module
USE Narran_Module
USE totj_Module
USE Testing_Module
USE FileUnits_Module
!
!            P U R P O S E   O F    S U B R O U T I N E
! This routine determines the overlap matrix elements between surface
! functions at different rho values.
!            I N P U T    A R G U M E N T S
!   rholast  previous hyperradius.
!   rho      current hyperradius
!   nang     number of quadrature angles.
!   nbasis   number of raw primitive basis functions.
!   nbasiss  number of symmetry-adapted primitive basis functions.
!   narran   number of arrangement channels.
!   xchanl    an array giving the channel number of each basis
!            function.
!   nangch   an array containing the starting position of the
!            angles for each of the arrangement channels.
!            nangch(narran+1) is equal to nang+1.
!   Phi_FEM  primitive basis of previous rho at present quadrature pts.
!   Phi_ABM   primitive basis of present rho at present quadrature pts.
!   vect     coeffs of present surface functions.
!            O U T P U T    A R G U M E N T S
!   s        overlap matrix elements.

IMPLICIT NONE
CHARACTER(LEN=200) :: filename, prefix
!            I N T E G E R S

INTEGER :: nang, ibasis, jbasis, nbasis, iang, iangst1, iangend1, xchanl, nangch, jpar, neigmin
INTEGER :: lnblnk, nbasiss, jrot, kbasis, nmode, irhfem, nprimabm, nprimfem, nabm, nfem, ndimen
!            R E A L S

REAL(Kind=WP_Kind) :: overlap, Phi_FEM, Phi_ABM, s, rholast, rho, vect,amatrx, xnorm, ThetAPH, chiaph, weight
!            D I M E N S I O N S
DIMENSION Phi_FEM (nang, nbasis), s (nbasiss, nbasiss), xchanl (nbasis), nangch (narran + 1)
DIMENSION Phi_ABM (nang, nbasis), vect (nbasiss, nbasiss), amatrx (nbasiss, nbasiss), jrot (nbasis)
DIMENSION ThetAPH (nang), chiaph (nang), weight (nang)

REAL(Kind=WP_Kind) :: diagonal (1000)
!            E X T E R N A L S

EXTERNAL popt, mxoutf, sgemm_junk, rdwrsurb
!-----------------------------------------------------------------------
! Determine printing options.
!-----------------------------------------------------------------------
LOGICAL :: little, medium, full, there
INTEGER :: ithcall, ithsub
DATA ithcall / 0 /, ithsub / 0 /
DATA little / .false. /, medium / .false. /, full / .false. /

CALL popt ('ovrfem  ', little, medium, full, ithcall, ithsub)
ndimen = nbasiss

INQUIRE(file = '.FEM_to_ABM', exist = there)
IF(there)THEN
   filename = '.FEM_to_ABM'
ELSE
   CALL getenv ('APH3D_HOME', prefix)
   filename = prefix (:lnblnk (prefix) ) //'/default_data/.FEM_to_ABM '
ENDIF

OPEN(Unit=FEM_TO_ABM_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//filename, status = 'old')
READ(FEM_TO_ABM_Unit, testing)
CLOSE(unit=FEM_TO_ABM_Unit)

IF(test=='ABM')THEN
   WRITE(Msg_Unit, * ) 'Testing ABM Surface Functions'
   WRITE(Out_unit, * ) 'Testing ABM Surface Functions'
ELSEIF (test=='FEM')THEN
   WRITE(Msg_Unit, * ) 'Testing FEM Surface Functions'
   WRITE(Out_unit, * ) 'Testing FEM Surface Functions'
ELSE
   WRITE(Msg_Unit, * ) 'Calculating FEM_TO_ABM transformation'
   WRITE(Out_unit, * ) 'Calculating FEM_TO_ABM transformation'

ENDIF
WRITE(Out_unit, * ) 'nangch=', nangch, ' nang=', nang
CALL basfem (ThetAPH, chiaph, nang, Phi_FEM, nmode, rholast, irhfem)
IF(symmetry)THEN
   CALL vscale (nang * ndimen, Phi_FEM, 1, One/SQRT(two) )
ENDIF
!-------------------------------------------------------------------
! Only used to isolate inaccuracies in ABM or FEM basis.
!-------------------------------------------------------------------
IF(test=='ABM')THEN
   WRITE(Msg_Unit, * ) 'Both Phi_FEM and Phi_ABM are ABM s.f.'
   WRITE(Out_unit, * ) 'Both Phi_FEM and Phi_ABM are ABM s.f.'
   CALL vcopy (nang * ndimen, Phi_FEM, Phi_ABM)
   nabm = neigmin
   nfem = neigmin
   nprimfem = nbasiss
   nprimabm = nbasiss
ELSEIF (test=='FEM')THEN
   WRITE(Msg_Unit, * ) 'Both Phi_FEM and Phi_ABM are FEM s.f.'
   WRITE(Out_unit, * ) 'Both Phi_FEM and Phi_ABM are FEM s.f.'
   CALL vcopy (nang * ndimen, Phi_ABM, Phi_FEM)
   nabm = nmode
   nfem = nmode
   nprimfem = nmode
   nprimabm = nmode
ELSE
   nabm = neigmin
   nfem = nmode
   nprimfem = nmode
   nprimabm = nbasiss

ENDIF
WRITE(Out_unit, * ) 'Number of ABM Surface Functions =', nabm
WRITE(Out_unit, * ) 'Number of FEM Surface Functions =', nfem
WRITE(Out_unit, * ) 'Number of Primitive ABM Basis Functions=', nprimabm

WRITE(Out_unit, * ) 'Number of Primitive FEM Basis Functions=', nprimfem
!
! Multiply the ABM and FEM surface functions by the weight.
!
CALL weight_set (nang, nprimfem, Phi_FEM, weight)

CALL weight_set (nang, nprimabm, Phi_ABM, weight)
! This will print weighted surface functions.
IF(medium)THEN
   WRITE(Out_unit, * ) 'chanl,iang,Phi_FEM,Phi_ABM'
   iangst1 = nangch (xchanl (1) )
   iangend1 = nangch (xchanl (1) + 1) - 1
   WRITE(Out_unit, * ) 'iangst1,iangend1', iangst1, iangend1
   DO 875 iang = iangst1, iangend1
      IF(abs (Phi_FEM (iang, 1) ) >0.05.or.abs (Phi_ABM (iang, 1) ) >0.05)THEN
         WRITE(Out_unit, * ) xchanl (1), iang, Phi_FEM (iang, 1), Phi_ABM (iang, 1)
      ENDIF
  875    ENDDO



ENDIF
!
!  --------------------------------------------------------------
!  evaluate overlaps of primitive bases of previous and present rho.
!  --------------------------------------------------------------
WRITE(Msg_Unit, * ) 'calculate overlap integrals'
WRITE(Out_unit, * ) 'calculate overlap integrals'
WRITE(Out_unit, * ) 'symmetry=', symmetry
IF(symmetry)THEN
   IF(jeven)THEN
      jpar = 0
   ELSE
      jpar = 1
   ENDIF
ENDIF
DO ibasis = 1, nprimabm
   iangst1 = nangch (xchanl (ibasis) )
   iangend1 = nangch (xchanl (ibasis) + 1) - 1
   WRITE(Out_unit, * ) 'ibasis=', ibasis, iangst1, iangend1, xchanl ( ibasis)
   DO jbasis = 1, nprimfem
      overlap = 0.0d0
      DO iang = iangst1, iangend1
         IF(symmetry.AND.xchanl (ibasis) ==2)THEN
            kbasis = ibasis + nbasis - nbasiss
            IF(mod (jrot (ibasis), 2) ==jpar)THEN
               overlap = overlap + Phi_FEM (iang, jbasis) * (Phi_ABM (iang, ibasis) + Phi_ABM (iang, kbasis) ) * sqrt (two)
            ELSE
               overlap = overlap + Phi_FEM (iang, jbasis) * (Phi_ABM (iang, ibasis) - Phi_ABM (iang, kbasis) ) * sqrt (two)
            ENDIF
         ELSE
            overlap = overlap + Phi_FEM (iang, jbasis) * Phi_ABM (iang, ibasis)
         ENDIF
      ENDDO
      amatrx (jbasis, ibasis) = overlap
   ENDDO
ENDDO
WRITE(Msg_Unit, * ) 'Overlap integrals have been calculated'

WRITE(Out_unit, * ) 'Overlap integrals have been calculated'
!--------------------------------------------------------------------
! Transform primitive overlaps integrals to actual overlaps.
!--------------------------------------------------------------------
WRITE(Out_unit, * ) 'Overlaps of Primitive Basis'
WRITE(Out_unit,  * ) 'Matrix dimensions:', nprimfem, '  by', nprimabm

CALL mxoutf (amatrx, nprimfem, nprimabm, ndimen, ndimen)
IF(test/='FEM')THEN
   WRITE(Out_unit, * ) 'Coefficients of the primative basis'
   CALL mxoutf (vect, nprimabm, nabm, ndimen, ndimen)

ENDIF
IF(test/='FEM')THEN
WRITE(Msg_Unit,  * ) 'Transform primitive overlaps to actual overlaps'
WRITE(Out_unit,  * ) 'Transform primitive overlaps to actual overlaps'
   CALL sgemm_junk (nprimfem, nprimabm, nabm, amatrx, ndimen, vect, ndimen, s, ndimen, 0, 1)
ELSEIF (test=='FEM')THEN
   CALL vcopy (ndimen * ndimen, s, amatrx)

ENDIF
IF(test=='ABM')THEN
   CALL sgemm_junk (nabm, nprimabm, nabm, vect, ndimen, s, ndimen, amatrx, ndimen, 4, 1)
   CALL vcopy (ndimen * ndimen, s, amatrx)

ENDIF
WRITE(Out_unit, * ) 'Overlaps Actual basis'
WRITE(Out_unit,  * ) 'Matrix dimensions:', nfem, '  by', nabm

CALL mxoutf (s, nfem, nabm, ndimen, ndimen)
!  --------------------------------------------------------------
!  WRITE overlap matrix on unit Ovrlp_FEM_Unit 'Ovrlp'
!  --------------------------------------------------------------
WRITE(Ovrlp_FEM_Unit) 0, rholast, rho
WRITE(Ovrlp_FEM_Unit) nabm, rholast, rho
DO jbasis = 1, nabm
   WRITE(Ovrlp_FEM_Unit) (s (ibasis, jbasis), ibasis = 1, nfem)
ENDDO
    
WRITE(Out_unit,'(A)') 'Diagonal Elements of Overlap basis'
DO ibasis = 1, min (nfem, nabm)
   diagonal (ibasis) = s (ibasis, ibasis)
ENDDO

WRITE(Out_unit,'(10ES11.4)') (diagonal (ibasis), ibasis = 1, min (nfem, nabm) )
!-----------------------------------------------------------------------
! IF desired WRITE out the overlap matrix.
!-----------------------------------------------------------------------
WRITE(Out_unit, * ) 'routine ovrbas'
WRITE(Out_unit, * ) 'nfem = ', nfem

WRITE(Out_unit, * ) 'nabm = ', nabm
WRITE(Out_unit, * ) 'Part of the overlap matrix at rholast, rho = ', rholast, rho
WRITE(Out_unit,  * ) 'Matrix dimensions:', nfem, '  by', nabm
CALL mxoutf (s, nfem, nabm, ndimen, ndimen)
!--------------------------------------------------------------------
!  construct s*s(trans) to test completeness of new basis.
!--------------------------------------------------------------------
WRITE(Out_unit,*) 'new basis completeness test.  rholast, rho= ', rholast, rho
CALL vsets (ndimen * ndimen, amatrx, 1, Zero)
CALL sgemm_junk (nfem, nabm, nfem, s, ndimen, s, ndimen, amatrx, ndimen, 2, 1)
CALL mxoutf (amatrx, nfem, nfem, ndimen, ndimen)
WRITE(Out_unit, * ) 'Diagonal Contributions'
DO ibasis = 1, nfem
   diagonal (ibasis) = amatrx (ibasis, ibasis)
ENDDO
WRITE(Out_unit, * ) (diagonal (ibasis), ibasis = 1, nfem)
!--------------------------------------------------------------------
!  construct s(trans)*s to test completeness of old basis.
!--------------------------------------------------------------------
WRITE(Out_unit,*) 'old basis completeness test.  rholast, rho= ', rholast, rho
CALL vsets (ndimen * ndimen, amatrx, 1, Zero)
CALL sgemm_junk (nfem, nabm, nabm, s, ndimen, s, ndimen, amatrx, ndimen, 4, 1)
CALL mxoutf (amatrx, nfem, nabm, ndimen, ndimen)
WRITE(Out_unit, * ) 'Diagonal Contributions'
DO ibasis = 1, nabm
   diagonal (ibasis) = amatrx (ibasis, ibasis)
ENDDO

WRITE(Out_unit, * ) (diagonal (ibasis), ibasis = 1, nabm)
IF(test=='ABM')THEN
   WRITE(Out_unit, * ) 'END of ABM test'
   STOP 'ABM test'
ELSEIF (test=='FEM')THEN
   WRITE(Out_unit, * ) 'END of FEM test'
   STOP 'FEM test'
ELSE
   WRITE(Out_unit, * ) 'Calculated Matrix between FEM and ABM ', 'Surface Functions'

ENDIF
RETURN
ENDSUBROUTINE ovrfem
