SUBROUTINE corgen(id, nelem, numnp, ntheta, nchi, Virtualx, idt, corner)
USE FileUnits_Module
!
!     $RCSfile: corgen.f,v $ $Revision: 1.14 $
!     $Date: 89/10/18 14:17:57 $
!     $State: Stable $
!


!     I N P U T    A R G U M E N T S

!     id(numnp)   id=-nce if the is constrained and positive otherwise
!     nelem
!     numnp    Number of nodal poin
!     ntheta
!     nchi
!     Virtualx
!     idt      Identifies each node as free, fixed, constrained, unique, nonexistant
!     corner


IMPLICIT NONE

!     L O G I C A L S

LOGICAL little, medium, full

!     I N T E G E R S
INTEGER id, nelem, numnp, corner, Virtualx, inode, iel, k, ntheta, nchi, itheta, jchi, imap, nel2, idt, ithcall, ithsub


!     D I M E N S I O N S

DIMENSION corner(1)
DIMENSION id(9,nelem), Virtualx(ntheta,nchi)
DIMENSION imap(9), idt(ntheta,nchi)

!     E X T E R N A L S
EXTERNAL popt
DATA imap/7,1,3,9,4,2,6,8,5/
DATA ithcall/0/, ithsub/0/
DATA little/.false./, medium/.false./, full/.false./
CALL popt ('corgen  ', little, medium, full, ithcall, ithsub)
DO jchi=1,nchi
   DO itheta=1,ntheta
      IF(idt(itheta,jchi)==-3)Virtualx(itheta,jchi)=0
   ENDDO
ENDDO
nel2=nelem/2
DO inode=1,numnp
   corner(inode)=0
ENDDO
DO iel=1,nel2
   DO k=1,4
      inode=id(imap(k),iel)
      corner(inode)=corner(inode)+2
   ENDDO
ENDDO
DO itheta=1,ntheta
   inode=Virtualx(itheta,1)
   IF(inode/=0)THEN
      corner(inode)=corner(inode)+4
   ENDIF
   inode=Virtualx(itheta,(nchi+1)/2)
   IF(inode/=0)THEN
      IF(corner(inode)==4.AND.corner(Virtualx(itheta,1))==8) corner(inode)=8
   ENDIF
ENDDO
DO jchi=1,(nchi+1)/2
   inode=Virtualx(1,jchi)
   IF(inode/=0)THEN
      corner(inode)=corner(inode)+4
   ENDIF
   inode=Virtualx(ntheta,jchi)
   IF(inode/=0)THEN
      corner(inode)=corner(inode)+4
   ENDIF
ENDDO
inode=Virtualx(1,1)
corner(inode)=8
inode=Virtualx(1,(nchi+1)/2)
corner(inode)=8
inode=Virtualx(ntheta,1)
corner(inode)=8
inode=Virtualx(ntheta,(nchi+1)/2)
corner(inode)=8
DO inode=1,numnp
   IF(corner(inode)==8)THEN
      corner(inode)=1
   ELSE
      corner(inode)=0
   ENDIF
ENDDO
IF(medium)THEN
   DO jchi=1,(nchi+1)/2
      DO itheta=1,ntheta
         inode=Virtualx(itheta,jchi)
         IF(corner(inode)==1)WRITE(Out_unit,*)'corner at: itheta,jchi,','node: ',itheta,jchi,inode
      ENDDO
   ENDDO
ENDIF
RETURN
END
