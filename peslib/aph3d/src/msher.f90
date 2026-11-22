SUBROUTINE msher (rhoval, mxelmnt, mxnode, mxtheta, mxchi, lnod, idt, scr1, scr2, iscr1, iscr2, virtualx,    &
           nid, beta, nodj, id, th, ch, idteq, idiv, lcor, mapr, ido, idivo, lcoro, mapnod,                  &
           sinlist, philist, iwork33, ithrho, numnp,  ndisce, nidm, neq, ndiv, nel, ntheta, nchi)
USE Numeric_Kinds_Module
USE Narran_Module
USE Converge_Module
USE mdata_Module
USE rhosur_Module
USE totj_Module
USE Masses_Module
USE convrsns_Module
USE elemnt_Module
USE FileUnits_Module
USE Numbers_Module
USE Point_RmsPhi_Module
!
!     $RCSfile: msher.f,v $ $Revision: 1.14 $
!     $Date: 89/10/18 14:18:10 $
!     $State: Stable $
!


!     I N P U T    A R G U M E N T S
!     rhoval
!     mxelmnt
!     mxnode
!     mxtheta
!     mxchi
!     lnod
!     idt      Identifies each node as free, fixed, constrained, unique, nonexistant
!     scr1
!     scr2
!     iscr1
!     iscr2
!     virtualx
!     nid(NDISCE) Array with the number of terms in the constraint equations
!     beta(ndisce) Constraint Coefficients
!     nodj
!     id(numnp)   id=-nce if the is constrained and positive otherwise
!     th
!     ch
!     idteq
!     idiv
!     lcor
!     mapr
!     ido
!     idivo
!     lcoro
!     mapnod
!     sinlist
!     philist
!     phirms
!     iwork33
!     ithrho
!     numnp    Number of nodal poin
!     ndisce    number of displacement constraints.
!     nidm
!     neq
!     ndiv
!     nel
!     ntheta
!     nchi


IMPLICIT NONE
LOGICAL little, medium, full, mapem
INTEGER i, iel, ithrho, jmax, nchi, ndisce, ndiv, nelem
INTEGER nident, nidm, ntheta, numnod, numnp, neq, ithcall, ithsub
INTEGER mxelmnt, mxnode, mxtheta, mxchi, idivo(mxelmnt), virtualx(mxtheta,mxchi)
INTEGER nel, nfixed, lnod(2*mxelmnt), idt(mxtheta*mxchi)
INTEGER iscr1(mxtheta*(mxchi/2+1)), iscr2(mxelmnt), nid(2*mxnode/3), nodj(nidm+1,2*mxnode/3+1)
INTEGER id(9,mxelmnt), idteq(mxnode), idiv(mxelmnt), lcor(2*mxelmnt), mapr(mxelmnt), ido(9*mxelmnt)
INTEGER lcoro(2*mxelmnt), mapnod(mxnode), iwork33(2*mxelmnt)
REAL(Kind=WP_Kind) chimax, chisym, scr1(mxnode), scr2(mxnode), rhoval
REAL(Kind=WP_Kind) sinlist(mxelmnt), philist(mxelmnt)       !, phirms(ntheta,nchi)
REAL(Kind=WP_Kind) beta(nidm,2*mxnode/3+1), th(9,mxelmnt), ch(9,mxelmnt)
!-----------------------------------------------------------------------
!     default DATA for 3d reactive scattering
!-----------------------------------------------------------------------
DATA mapem /.false./
DATA ithcall/0/, ithsub/0/
DATA little/.false./, medium/.false./, full/.false./
CALL popt ('msher   ', little, medium, full, ithcall, ithsub)
!-----------------------------------------------------------------------
!     NSYM = 1  Heteronuclear diatomic in channel 1.
!     NSYM = 2  Homonuclear diatomic in channel.
!-----------------------------------------------------------------------
IF(nsymc(1))THEN
   symmetry = .true.
ELSE
   symmetry = .false.
ENDIF
rhosurf=rhoval
IF(ndiv>10)THEN
   WRITE(Out_unit,*) '***error***:ndiv>10:',ndiv
   STOP 'msher'
ENDIF
IF(symmetry.AND.(mass(2)/=mass(3)))THEN
   WRITE(Out_unit,*)'***error***:symmetry.AND.(mass(2)/=mass(3)):', symmetry,mass(2),mass(3)
   STOP 'msher'
ENDIF
!-----------------------------------------------------------------------
!     find the largest chi and the symmetry chi
!-----------------------------------------------------------------------
IF(.NOT.symmetry) nident=2
IF(symmetry) nident=4
chisym=eight*atan(one)/nident
chimax=2*chisym
WRITE(Out_unit,*) 'starting rhosurf=',rhosurf,' parity=',parity,' jtot=',jtot,' mega=',mega
!-----------------------------------------------------------------------
!     determine the number of rows and columns.
!-----------------------------------------------------------------------
!     ntheta=nlt*(2**ndiv)*2+1
jmax=0
DO i=1,3
   jmax=jmax+nlc(i)*(2**ndiv)*2+1
   IF(i/=1)jmax=jmax-1
ENDDO
IF(symmetry)jmax=(jmax+1)/2
!     nchi=2*jmax-1
IF(ntheta*nchi>mxtheta*mxchi)THEN
   WRITE(Out_unit,*)'error ntheta,nchi,mxelmnt = ',ntheta,nchi,mxelmnt
   STOP 'msher'
ENDIF
!-----------------------------------------------------------------------
!     determine the mesh.
!-----------------------------------------------------------------------
CALL mesher(nlt, nlc, ndiv, id, th, ch, mxelmnt, numnp,  nelem, chisym, mapem, nodmax, idiv, numnod, lnod,               &
           IntTheta, IntChi, lcor, scr1, scr2, iscr1, iscr2, virtualx, mxtheta, mxchi, mapr, ido, idivo, lcoro, mapnod, sinlist,    &
           philist, phirms, mxnode, scr1, scr2, iwork33, ithrho, ntheta, nchi)
!-----------------------------------------------------------------------
!     Set up the constraint equations.
!-----------------------------------------------------------------------
CALL cnstrnt(ntheta,nchi,idt, lnod, parity, jeven, symmetry, nid, nodj, beta, nidm, ndisce, numnp, mxelmnt, jmax, nelem, id,  &
           scr1, scr2, virtualx, iscr1, nfixed, mega, IntTheta, IntChi)
!-----------------------------------------------------------------------
!     WRITE DATA for the surface function calculation.
!-----------------------------------------------------------------------
nmode=nfreq(0)
WRITE(Msher_Bin_Unit)  nelem, numnp, ndisce, nidm, nfreq, nq, nmode, rtol
!-----------------------------------------------------------------------
!     WRITE the constraint equations.
!-----------------------------------------------------------------------
CALL nWRITE(ntheta, nchi, idt, numnp, nodj, beta, nid,  nidm, ndisce, idteq, scr1, scr2, virtualx)
neq=numnp-ndisce-nfixed
!-----------------------------------------------------------------------
!     WRITE the nodal DATA.
!     The order written is prescribed by the surface function routines.
!-----------------------------------------------------------------------
DO iel=1,nelem
   IF(full)THEN
      WRITE(Out_unit,*) iel,id(1,iel),id(3,iel),id(9,iel),id(7,iel),id(2,iel),id(6,iel),id(8,iel),id(4,iel),id(5,iel)
   ENDIF
   WRITE(Msher_Bin_Unit) iel, id(1,iel), id(3,iel), id(9,iel), id(7,iel), id(2,iel), id(6,iel), id(8,iel), id(4,iel), id(5,iel)
ENDDO
nel=nelem

RETURN
END
