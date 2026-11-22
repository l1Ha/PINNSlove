SUBROUTINE hiphi(nnc, chii, dthe, dchi, jbox, ibox, iel,lcor, linkel, jreg, sinlist, philist, phirms, &
           mxnode, ntheta, nchi, mxelmnt, lowel, idiv)
USE Numeric_Kinds_Module
USE Numbers_Module
!
!     $RCSfile: hiphi.f,v $ $Revision: 1.14 $
!     $Date: 89/10/18 14:18:04 $
!     $State: Stable $
!

!     P U R P O S E   O F    S U B R O U T I N E

!     This routine replaces old routine lowpot.
!     This routine finds the maximum phirms for the 4 elements that have
!     just been created from one old element; this should be run just after
!     divide been run.


!     I N P U T    A R G U M E N T S
!
!     nnc
!     chii
!     dthe
!     dchi
!     jbox
!     ibox
!     iel
!     lcor
!     linkel
!     jreg
!     sinlist
!     philist
!     phirms
!     mxnode
!     ntheta
!     nchi
!     mxelmnt
!     lowel
!     idiv


IMPLICIT NONE
INTEGER :: i, ibox, icor, iel, j, jj, jbox, jcor, jel,  jreg, mxelmnt, mi, mj, nchi,  mxnode, ntheta
INTEGER :: idiv(mxelmnt), lcor(2,mxelmnt), linkel(mxelmnt), lowel(mxelmnt), nnc(3)
REAL(Kind=WP_Kind) :: dchi(3), chii(3), philist(mxelmnt), phirms(ntheta,nchi), sinlist(mxelmnt)
REAL(Kind=WP_Kind) :: chi, dthe, phim, theta, thetm

jel = iel
DO 175 jj=1,4
   icor = lcor(1,jel)
   jcor = lcor(2,jel)

   phim=zero
   DO 166 mj=0,2*jbox
      j=jcor+mj
      IF(jreg==1)THEN
         chi=(j-1)*dchi(jreg)
      ELSEIF (jreg==2)THEN
         chi=chii(2)+(j-nnc(1))*dchi(jreg)
      ELSE
         chi=chii(3)+(j-nnc(2)-nnc(1)+1)*dchi(jreg)
      ENDIF
      DO 167 mi=0,2*ibox
         i=icor+mi
         theta=(i-1)*dthe
         IF(phirms(i,j)>phim)THEN
            phim=phirms(i,j)
            thetm=theta
         ENDIF
 167        CONTINUE
 166     CONTINUE
   philist(jel) = phim
   sinlist(jel)=sin(thetm)
   IF(philist(jel)<=zero)philist(jel)=1.d-60
   IF(sinlist(jel)<=zero)sinlist(jel)=1.d-20
   jel = linkel(jel)
 175  CONTINUE

RETURN
END
