SUBROUTINE divide(mxelmnt, nelem, jbox, ibox, imax, jdiv, iel, icor, jcor, nelnew, numnod, idiv, lcor, lnod, linkel, numsto)
USE FileUnits_Module
!
!     $RCSfile: divide.f,v $ $Revision: 1.14 $
!     $Date: 89/10/18 14:17:59 $
!     $State: Stable $
!


!     I N P U T    A R G U M E N T S
!
!     mxelmnt
!     nelem
!     jbox
!     ibox
!     imax
!     jdiv
!     iel
!     icor
!     jcor
!     nelnew
!     numnod
!     idiv
!     lcor
!     lnod
!     linkel
!     numsto


IMPLICIT NONE

!     I N T E G E R S
INTEGER ibcor, ibox, icor, idiv, iel, imax, jbcor, jcor, jdiv, lcor, linkel, lnod, mxelmnt, mi, mj, nelem, nelnew, numnod, &
           numsto, jbox


!     D I M E N S I O N S
DIMENSION linkel(mxelmnt), lcor(2,mxelmnt), lnod(2*mxelmnt), idiv(mxelmnt)


!     E X T E R N A L S
EXTERNAL nodson



!-----------------------------------------------------------------------
!     routine to acutally divide an old element into 4 new elements and upda
!     the element link list
!
!     linkel, lnod, nelnew, numnod, idiv, lcor are changed on RETURN
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     divide the old element into 4 elements
!-----------------------------------------------------------------------
DO mj=0,1
   jbcor=jcor+mj*jbox
   DO mi=0,1
      ibcor=icor+mi*ibox
      CALL nodson (ibcor, jbcor, ibox, jbox, numnod, lnod, imax, mxelmnt, numsto)
      IF(mj>0.or.mi>0)THEN
         nelnew=nelnew+1
         idiv(nelem+nelnew)=jdiv
         lcor(1,nelem+nelnew)=ibcor
         lcor(2,nelem+nelnew)=jbcor
      ENDIF
    ENDDO
ENDDO
idiv(iel)=jdiv
!-----------------------------------------------------------------------
!     update the elements link list
!-----------------------------------------------------------------------
linkel(nelem+nelnew)=linkel(iel)
linkel(iel)=nelem+nelnew-2
linkel(nelem+nelnew-2)=nelem+nelnew-1
linkel(nelem+nelnew-1)=nelem+nelnew
RETURN
END
