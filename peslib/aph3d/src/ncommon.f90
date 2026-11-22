SUBROUTINE ncommon(ntheta, nchi, idt, nelem, id, IntTheta, IntChi)
USE FileUnits_Module
USE FEM_Module
!
!     $RCSfile: ncommon.f,v $ $Revision: 1.14 $
!     $Date: 89/10/18 14:18:12 $
!     $State: Stable $
!

!     P U R P O S E   O F    S U B R O U T I N E


!     I N P U T    A R G U M E N T S
!     ntheta
!     nchi
!     idt      Identifies each node as free, fixed, constrained, unique, nonexistant
!     nelem
!     id(numnp)   id=-nce if the is constrained and positive otherwise
!     IntTheta
!     IntChi
IMPLICIT NONE
INTEGER ntheta, nchi, jchi, itheta, node, ielem, nelem, k
INTEGER idt(ntheta,nchi), id(9,nelem), imap(9)
INTEGER IntTheta(9,nelem), IntChi(9,nelem)
!-----------------------------------------------------------------------
!     This routine finds nodes that are not in common with adjacent
!     elements and hence must be constrained.
!-----------------------------------------------------------------------
DATA imap/7,1,3,9,4,2,6,8,5/

DO ielem=1,nelem
   DO k=1,8
      node=id(imap(k),ielem)
      itheta=IntTheta(imap(k),ielem)
      jchi=IntChi(imap(k),ielem)
      IF(idt(itheta,jchi)>=0.AND.idt(itheta,jchi)/=1)THEN
         idt(itheta,jchi)=idt(itheta,jchi)+2
      ENDIF
   ENDDO
ENDDO
!-----------------------------------------------------------------------
!     The side nodes should be in common with a node of one other element.
!-----------------------------------------------------------------------
DO ielem=1,nelem
   DO 30 k=5,8
      itheta=IntTheta(imap(k),ielem)
      jchi=IntChi(imap(k),ielem)
      IF(idt(itheta,jchi)==2)THEN
         IF(itheta==1)GOTO 30
         IF(itheta==ntheta)GOTO 30
         IF(jchi==1)GOTO 30
         IF(jchi==nchi)GOTO 30
         idt(itheta,jchi)=UniqueNode
      ENDIF
 30      CONTINUE
ENDDO
!-----------------------------------------------------------------------
!     The corner nodes should be in common with a node from four other
!     elements.
!-----------------------------------------------------------------------
DO 4 ielem=1,nelem
   DO 40 k=1,4
      itheta=IntTheta(imap(k),ielem)
      jchi=IntChi(imap(k),ielem)
      IF(idt(itheta,jchi)>=0.AND.idt(itheta,jchi)/=1)THEN
         IF(itheta==1)GOTO 4
         IF(itheta==ntheta)GOTO 4
         IF(jchi==1)GOTO 4
         IF(jchi==nchi)GOTO 4
         IF(idt(itheta,jchi)==8)THEN
            GO TO 40
         ELSEIF(idt(itheta,jchi)==6)THEN
            GO TO 40
         ELSEIF(idt(itheta,jchi)==4)THEN
            idt(itheta,jchi)=UniqueNode
            GO TO 40
         ELSEIF(idt(itheta,jchi)==2)THEN
            idt(itheta,jchi)=UniqueNode
            GO TO 40
         ELSE
            WRITE(Out_unit,*)'***error*** corner node problem'
            WRITE(Out_unit,*)'itheta,jchi,ielem,idt = ',itheta, jchi,ielem,idt(itheta,jchi)
            STOP 'ncommon'
         ENDIF
      ENDIF
 40      CONTINUE
 4    CONTINUE
DO jchi=1,nchi
   DO itheta=1,ntheta
      IF(idt(itheta,jchi)>=0.AND.idt(itheta,jchi)/=1)THEN
         idt(itheta,jchi)=0
      ENDIF
   ENDDO
ENDDO
!     temporary reads
REWIND TempNode_Bin_Unit
DO ielem=1,nelem
   DO k=1,9
      WRITE(TempNode_Bin_Unit)id(k,ielem),IntTheta(k,ielem),IntChi(k,ielem)
   ENDDO
ENDDO
RETURN
END
