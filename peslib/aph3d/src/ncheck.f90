SUBROUTINE ncheck(ntheta, nchi, idt, numnp, nfree, ndisce, nfixed)
USE Numeric_Kinds_Module
USE FileUnits_Module
USE FEM_Module
USE Numbers_Module

!     I N P U T    A R G U M E N T S
!     ntheta
!     nchi
!     idt      Identifies each node as free, fixed, constrained, unique, nonexistant

!     O U T P U T    A R G U M E N T S
!     nfree     Number of free nodes
!     ndisce    Number of displacement constraints.
!     nfixed    Number of nodes fixed to be zero
!     numnp     Number of nodal points


IMPLICIT NONE
INTEGER, INTENT(IN)  :: ntheta, nchi
INTEGER, INTENT(IN)  :: idt(ntheta,nchi) 
INTEGER, INTENT(OUT) :: nfree, ndisce, nfixed, numnp
INTEGER              :: jchi, itheta, idtval, nonexist, nunique, ntotal

ntotal=ntheta*nchi
ndisce=0
nfree=0
nfixed=0
nonexist=0
nunique=0
DO jchi=1,nchi
   DO itheta=1,ntheta
      idtval=idt(itheta,jchi)
      IF(idtval==nonxst)THEN
         nonexist=nonexist+1
      ELSEIF(idtval==cnstrain)THEN
         ndisce=ndisce+1
      ELSEIF(idtval==free)THEN
         nfree=nfree+1
      ELSEIF(idtval==fixed)THEN
         nfixed=nfixed+1
      ELSEIF(idtval==UniqueNode)THEN
         nunique=nunique+1
      ELSE
         WRITE(Out_unit,*)'error:',itheta,jchi,idtval
         STOP 'ncheck'
      ENDIF
   ENDDO
ENDDO
numnp=nfree+ndisce+nfixed
WRITE(Out_unit,FMT=2)ntotal
WRITE(Out_unit,FMT=3)nonexist
WRITE(Out_unit,FMT=4)ndisce
WRITE(Out_unit,FMT=5)nfixed
WRITE(Out_unit,FMT=6)numnp
WRITE(Out_unit,FMT=7)nfree
WRITE(Out_unit,FMT=8)nunique
IF(ndisce>int(numnp/ThreeHalves))THEN
   WRITE(Out_unit,*)'too many constrained nodes, ndisce=',ndisce
   WRITE(Out_unit,*)'greater that ',int(numnp/ThreeHalves)
   STOP 'ncheck'
ENDIF
 2    FORMAT(1x,'Total number of mesh points=',i6)
 3    FORMAT(1x,'Number of points that are not used=',i6)
 4    FORMAT(1x,'Number of constraints=',i5)
 5    FORMAT(1x,'Number of nodal values fixed at zero=',i5)
 6    FORMAT(1x,'Number of nodal points=',i5)
 7    FORMAT(1x,'Number of free nodes = number of equations=',i5)
 8    FORMAT(1x,'Number of unique nodes=',i5)
RETURN
END
