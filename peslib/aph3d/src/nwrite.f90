SUBROUTINE nWRITE(ntheta, nchi, idt, numnp, nodj, beta, nid, nidm, ndisce, idteq, th, ch, Virtualx)
USE Numeric_Kinds_Module
USE FileUnits_Module
!
!     $RCSfile: nWRITE.f,v $ $Revision: 1.14 $
!     $Date: 89/10/18 14:18:13 $
!     $State: Stable $
!


!     I N P U T    A R G U M E N T S
!
!     ntheta
!     nchi
!     idt      Identifies each node as free, fixed, constrained, unique, nonexistant
!     numnp    Number of nodal poin
!     nodj
!     beta(ndisce) Constraint Coefficients
!     nid(NDISCE) Array with the number of terms in the constraint equations
!     nidm
!     ndisce    number of displacement constraints.
!     idteq
!     th
!     ch
!     Virtualx


IMPLICIT NONE

!     L O G I C A L S

LOGICAL little, medium, full

!     I N T E G E R S

INTEGER cnstrain, nonxst, ithcall, ithsub
INTEGER ntheta, nchi, idt, numnp, nodj, nid, nidm, ndisce, inode, itheta, jchi, nce, nc, Virtualx, idteq

!     R E A L S
REAL(Kind=WP_Kind) th, ch, beta

!     D I M E N S I O N S

DIMENSION idt(ntheta,nchi), idteq(numnp), Virtualx(ntheta,nchi)
DIMENSION nodj(nidm+1,ndisce), beta(nidm,ndisce), nid(ndisce)
DIMENSION th(numnp), ch(numnp)

!     E X T E R N A L S

EXTERNAL recon, popt



PARAMETER (cnstrain=-2)
PARAMETER (nonxst=-3)
DATA ithcall/0/, ithsub/0/
DATA little/.false./, medium/.false./, full/.false./
CALL popt ('nwrite  ', little, medium, full, ithcall, ithsub)
DO 10 jchi=1,nchi
   DO 1 itheta=1,ntheta
      IF(idt(itheta,jchi)/=nonxst)THEN
         inode=Virtualx(itheta,jchi)
         IF(inode>=1.AND.inode<=numnp)THEN
            idteq(inode)=idt(itheta,jchi)
         ELSE
            WRITE(Out_unit,*)'itheta,jchi,inode= ',itheta,jchi,inode
            STOP 'nWRITE'
         ENDIF
      ENDIF
 1       CONTINUE
 10   CONTINUE
CALL recon (idteq, nodj, beta, nid, nidm, ndisce)
CALL recon (idteq, nodj, beta, nid, nidm, ndisce)
CALL recon (idteq, nodj, beta, nid, nidm, ndisce)
DO inode=1,numnp
   WRITE(Msher_Bin_Unit)inode, idteq(inode), th(inode), ch(inode)
   IF(full)THEN
      WRITE(Out_unit,'("inode=",i5," idteq=",i5," th=",es15.7," ch=",es15.7)')inode,idteq(inode),th(inode),ch(inode)
   ENDIF
ENDDO
DO 2 nce=1,ndisce
   DO 7 nc=1,nid(nce)
      inode=nodj(nc,nce)
      IF(inode>numnp)THEN
         WRITE(Out_unit,*)'error: inode>numnp ',inode,numnp
         WRITE(Out_unit,*)nce,nc,nid(nce),nodj(nc,nce),nidm, nodj(nidm+1,nce),beta(nc,nce)
         STOP 'nWRITE'
      ENDIF
      IF(idteq(nodj(nc,nce))/=0)GOTO 4
      IF(nodj(nc,nce)<=0.or.nodj(nc,nce)>numnp)GOTO 4
    7    CONTINUE
   IF(nodj(nidm+1,nce)<=0.or.nodj(nidm+1,nce)>numnp)GOTO 4
   WRITE(Msher_Bin_Unit)nce,nodj(nidm+1,nce),nid(nce)
   WRITE(Msher_Bin_Unit)(nodj(nc,nce),beta(nc,nce),nc=1,nid(nce))
   IF(medium)THEN
      WRITE(Out_unit,'("nce=",i5," nodj=",i5," nid=",i5)')nce,nodj(nidm+1,nce),nid(nce)
      WRITE(Out_unit,'("nodj=",i5," beta=",es15.7)')(nodj(nc,nce),beta(nc,nce),nc=1,nid(nce))
   ENDIF
    2 CONTINUE
RETURN
!-----------------------------------------------------------------------
!     This is reached only IF an error was found.
!-----------------------------------------------------------------------
 4    CONTINUE
DO 3 nce=1,ndisce
   DO 30 nc=1,nid(nce)
      inode=nodj(nc,nce)
      IF(idteq(nodj(nc,nce))/=0.or. nodj(nc,nce)<=0.or.nodj(nc,nce)>numnp.or. &
                 nodj(nidm+1,nce)<=0.or.nodj(nidm+1,nce)>  numnp)THEN
         WRITE(Out_unit,*)'error'
         WRITE(Out_unit,*)'nodj(nidm+1,nce),nid(nce)'
         WRITE(Out_unit,*)nodj(nidm+1,nce),nid(nce)
         WRITE(Out_unit,*)'nce,nc,inode,idteq(inode),th(inode),', 'ch(inode)'
         WRITE(Out_unit,*)nce,nc,inode,idteq(inode),th(inode), ch(inode)
      ENDIF
 30      CONTINUE
 3    CONTINUE
STOP 'nWRITE'
END
