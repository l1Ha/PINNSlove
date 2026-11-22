SUBROUTINE recon(idteq, nodj, beta, nid, nidm, ndisce)
USE Numeric_Kinds_Module
USE FileUnits_Module
!
!     $RCSfile: recon.f,v $ $Revision: 1.14 $
!     $Date: 89/10/18 14:18:15 $
!     $State: Stable $


IMPLICIT NONE

!     L O G I C A L S

LOGICAL little, medium, full

!     I N T E G E R S
INTEGER idteq, nodj, nid, nidm, ndisce, nce, nc, inode, ncep, inodep, ncp, ithcall, ithsub


!     R E A L S
REAL(Kind=WP_Kind) beta

!     D I M E N S I O N S

DIMENSION nodj(nidm+1,ndisce), beta(nidm,ndisce),  nid(ndisce), idteq(1)

!     E X T E R N A L S

EXTERNAL popt
DATA ithcall/0/, ithsub/0/
DATA little/.false./, medium/.false./, full/.false./
CALL popt ('recon   ', little, medium, full, ithcall, ithsub)
IF(little)THEN
   WRITE(Out_unit,*)'Inside recon: ndisce',ndisce
   WRITE(Out_unit,*)'Inside recon: nid',nid
   WRITE(Out_unit,*)'Inside recon: nodj',nodj
ENDIF
DO 20 nce=1,ndisce
   DO 10 nc=1,nid(nce)
      inode=nodj(nc,nce)
      IF(idteq(inode)==0)GOTO 10
      IF(medium)THEN
         WRITE(Out_unit,*)'node= ',nodj(nidm+1,nce), ' is constrained to ',inode
      ENDIF
      DO 5 ncep=1,ndisce
         inodep=nodj(nidm+1,ncep)
         IF(inode/=inodep)GOTO 5
         IF(medium)THEN
            WRITE(Out_unit,*)'but node = ',inodep, ' is constrained to node'
            WRITE(Out_unit,*)(nodj(ncp,ncep),ncp=1,nid(ncep))
         ENDIF
         DO 4 ncp=2,nid(ncep)
            nodj(nid(nce)+ncp-1,nce)=nodj(ncp,ncep)
            beta(nid(nce)+ncp-1,nce)=beta(nc,nce) *beta(ncp,ncep)
 4             CONTINUE
         nid(nce)=nid(nce)+nid(ncep)-1
         nodj(nc,nce)=nodj(1,ncep)
         beta(nc,nce)=beta(nc,nce)*beta(1,ncep)
         IF(medium)THEN
            WRITE(Out_unit,*)'nodj(nidm+1,nce) = ',nodj(nidm+1,nce)
            WRITE(Out_unit,*)'is now constrained to'
            WRITE(Out_unit,*)'nodj=',(nodj(ncp,nce),ncp=1,nid(nce))
            WRITE(Out_unit,*)'beta=',(beta(ncp,nce),ncp=1,nid(nce))
         ENDIF
         GOTO 10
 5          CONTINUE
 10      CONTINUE
 20   CONTINUE
RETURN
END
