SUBROUTINE rsffem (th, ch, nel, id, f, beta, idi, nid, nidm, freq, xnode, idni, nmode, numnp, ndisce, neq, symmetry)
USE Numeric_Kinds_Module
USE FileUnits_Module
!-----------------------------------------------------------------------
! Reads in Finite Element DATA
!-----------------------------------------------------------------------

IMPLICIT NONE
INTEGER :: id, nid, idi, idni, ndisce, nidm, numnp, nmode, neq, i, j, k, iel, imode, inode, nel, nelm, ii, nsyml
REAL(Kind=WP_Kind) :: th, ch, f, beta, freq, xnode, fac1
LOGICAL :: symmetry
DIMENSION id (9, nel), th (9, nel), ch (9, nel)

DIMENSION f (numnp, nmode), nid (ndisce), beta (nidm, ndisce), idi (nidm, ndisce), freq (nmode), xnode (numnp), idni (numnp)
WRITE(Msg_Unit, * ) 'Reading in FEM DATA'
READ(Nmodes_FEM_Bin_Unit, END = 90) (freq (i), i = 1, nmode)
WRITE(Out_Unit, * ) 'Finite Element Eigenenergies'
WRITE(Out_Unit, * ) freq
READ(Nmodes_FEM_Bin_Unit) (nid (i), i = 1, ndisce)
READ(Nmodes_FEM_Bin_Unit) ( (idi (i, j), i = 1, nidm), j = 1, ndisce)
READ(Nmodes_FEM_Bin_Unit) ( (beta (i, j), i = 1, nidm), j = 1, ndisce)
READ(Nmodes_FEM_Bin_Unit) (idni (i), i = 1, numnp)
!-----------------------------------------------------------------------
! READ in element definition DATA
!-----------------------------------------------------------------------
READ(ElmDef_Unit, END = 96) nelm, nsyml
READ(ElmDef_Unit, END = 95) ( (id (k, iel), k = 1, 9), iel = 1, nelm)
READ(ElmDef_Unit, END = 95) ( (th (k, iel), k = 1, 9), iel = 1, nelm)
READ(ElmDef_Unit, END = 95) ( (ch (k, iel), k = 1, 9), iel = 1, nelm)
!  ---------------------------------------------------------------------
!  READ surface functions
!  ---------------------------------------------------------------------
DO imode = 1, nmode
   READ(Nmodes_FEM_Bin_Unit, END = 90) ii, (f (i, imode), i = 1, neq)
!-----------------------------------------------------------------------
! reconstruct the surface function for the whole space
! condis sorts the surface fn. by eqn. numbers, so we must resort
! it by global node numbers.
!-----------------------------------------------------------------------
!     WRITE(Out_Unit,*) numnp,nidm,ndisce,neq
   CALL condis (nid, idi, beta, f (1, imode), nidm, ndisce, neq)
   DO inode = 1, numnp
      IF(idni (inode) <0)THEN
         xnode (inode) = f (neq - idni (inode), imode)
      ELSEIF (idni (inode) >0)THEN
         xnode (inode) = f (idni (inode), imode)
      ENDIF
   ENDDO
   DO inode = 1, numnp
      f (inode, imode) = xnode (inode)
   ENDDO
ENDDO
!-----------------------------------------------------------------------
! determine normalization factor and normalize surface function.
!-----------------------------------------------------------------------
fac1 = 1.0
IF(.NOT.symmetry)THEN
   fac1 = 1.0
!        fac1=1.0/sqrt(2.0)
ELSEIF (symmetry)THEN
   fac1 = 1.0 / sqrt (2.0)
ENDIF
DO inode = 1, numnp
   DO j = 1, nmode
      f (inode, j) = f (inode, j) * fac1
   ENDDO
ENDDO
!-----------------------------------------------------------------------
! f now contains a surface-function which is orthonormal.
!-----------------------------------------------------------------------
WRITE(Msg_Unit, * ) 'Finished reading FEM DATA'
RETURN
! IF an END of file is found, come here.
   90 WRITE(Out_Unit, * ) 'END of surface function file'
STOP 'rsffem'
   95 WRITE(Out_Unit, * ) 'END of elmdef file'
STOP 'rsffem'
   96 WRITE(Out_Unit, * ) 'END of elmdef file, first'
STOP 'rsffem'
!
   80 FORMAT(/,'there are ',i5,'  normal modes')
  140 FORMAT(/,'jcap=',i5,'  lambda=',i5,'  rho=',1pe21.14)
ENDSUBROUTINE rsffem
