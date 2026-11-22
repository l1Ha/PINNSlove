SUBROUTINE cnstrnt(ntheta, nchi, idt, lnod, parity, jeven, symmetry, nid, nodj, beta, nidm, ndisce, &
           numnp, mxelmnt, jmax, nelem, id, scr1, scr2, Virtualx, corner, nfixed, mega, IntTheta, IntChi)
USE Numeric_Kinds_Module
USE FileUnits_Module
USE FEM_Module
!
!     $RCSfile: cnstrnt.f,v $ $Revision: 1.14 $
!     $Date: 89/10/18 14:17:56 $
!     $State: Stable $
!


!     I N P U T    A R G U M E N T S
!
!
!     ntheta
!     nchi
!     idt
!     lnod
!     parity
!     jeven
!     symmetry
!     nid(NDISCE) Array with the number of terms in the constraint equations
!     nodj
!     beta(ndisce) Constraint Coefficients
!     nidm
!     ndisce    number of displacement constraints.
!     numnp    Number of nodal poin
!     mxelmnt
!     jmax
!     nelem
!     id(numnp)   id=-nce if the is constrained and positive otherwise
!     IntTheta
!     IntChi
!     scr1
!     scr2
!     Virtualx
!     corner
!     nfixed
!     mega

IMPLICIT NONE
LOGICAL jeven, symmetry
INTEGER idt, ntheta, nchi, jchi, itheta, nodj, nidm, ndisce, mega
INTEGER numnp, lnod, jmax, mxelmnt, nid, nelem, id, nfree, nfixed, parity, node, corner, Virtualx
INTEGER IntTheta(9,nelem), IntChi(9,nelem)          
REAL(Kind=WP_Kind) beta, scr1, scr2
DIMENSION idt(ntheta,nchi), lnod(numnp), corner(numnp), nid(numnp), nodj(nidm+1,2*numnp/3+1), beta(nidm,2*numnp/3+1), &
          id(numnp), Virtualx(ntheta,nchi)
DIMENSION scr1(numnp), scr2(numnp)
EXTERNAL  bitint, symnod, znodes, ncommon, eqncn2, ncheck, corgen, fixin

!-----------------------------------------------------------------------
!     Initialize idt array to nonxst
!     idt(itheta,jchi)=nonxst                implies a nonexistant node
!-----------------------------------------------------------------------
DO jchi=1,nchi
   DO itheta=1,ntheta
      idt(itheta,jchi)=nonxst
   ENDDO
ENDDO
!-----------------------------------------------------------------------
!     Determine free and nonexistant nodes
!     idt(itheta,jchi)= free                free node.
!     idt(itheta,jchi)= nonxst         implies a nonexistant node.
!-----------------------------------------------------------------------
CALL bitint(ntheta, nchi, idt, lnod, symmetry, jmax, 2*mxelmnt)
CALL corgen(id, nelem, numnp, ntheta, nchi, Virtualx, idt,  corner)
!-----------------------------------------------------------------------
! Count the total number of nodes
!-----------------------------------------------------------------------
node=0
DO jchi=1,nchi
   DO itheta=1,ntheta
      IF(idt(itheta,jchi)/=nonxst)node=node+1
   ENDDO
ENDDO
IF(node/=numnp)THEN
   WRITE(Msg_unit,*)'error: node,numnp',node,numnp
   WRITE(Out_unit,*)'error: node,numnp',node,numnp
   STOP 'cnstrnt'
ENDIF
!-----------------------------------------------------------------------
!     Apply symmetry constraints
!     idt(itheta,jchi)= cnstrain       node constrained to another node
!-----------------------------------------------------------------------
CALL symnod(ntheta, nchi, idt, nid, nodj, beta, nidm, ndisce, parity, jeven, symmetry, jmax, Virtualx, mega)
!-----------------------------------------------------------------------
!     Fix some node to be zero.
!     idt(itheta,jchi)= fixed          nodal value is zero.
!-----------------------------------------------------------------------
CALL znodes(ntheta, nchi, idt, parity, jeven, symmetry, mega)
!     CALL fixin(idt, ntheta, nchi, thetaval, chivals)
!-----------------------------------------------------------------------
!     Look for nodes that are not common to an adjacent element.
!     idt(itheta,jchi)= UniqueNode   node that is not common to an adjacent node.
!-----------------------------------------------------------------------
CALL ncommon(ntheta, nchi, idt, nelem, id, IntTheta, IntChi)
!-----------------------------------------------------------------------
!     Determine constraint equations.
!-----------------------------------------------------------------------
CALL eqncn2(ntheta, nchi, idt, nid, nodj, beta, nidm, ndisce, Virtualx, nelem, IntTheta, IntChi)
!-----------------------------------------------------------------------
!     Check the constraint conditions
!
!     idt(itheta,jchi)= fixed          nodal value is zero.
!     idt(itheta,jchi)= free           free node.
!     idt(itheta,jchi)= cnstrain       node constrained to another nod
!     idt(itheta,jchi)= nonxst         implies a nonexistant node.
!     idt(itheta,jchi)= UniqueNode         node that is not common to an
!     adjacent node.
!-----------------------------------------------------------------------
CALL ncheck(ntheta, nchi, idt, numnp, nfree, ndisce, nfixed)
RETURN
END
