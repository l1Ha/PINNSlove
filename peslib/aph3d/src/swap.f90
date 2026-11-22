SUBROUTINE swap (nelem, maxelp, lcor, idiv, id)
!
!     $RCSfile: swap.f,v $ $Revision: 1.14 $
!     $Date: 89/10/18 14:18:15 $
!     $State: Stable $
!

!     P U R P O S E   O F    S U B R O U T I N E
!-----------------------------------------------------------------------
!     this routine swaps two elements in lcor, idiv, & id
!     so that the point (imax,jmax) is in the element with the
!     largest element number
!-----------------------------------------------------------------------

IMPLICIT NONE

!     I N T E G E R S

INTEGER i, id, id1, id2, idiv, k, lcor, maxelp, nelem

!     D I M E N S I O N S

DIMENSION lcor(2,nelem), idiv(nelem), id(9,nelem)


DO 10 i=1,2
   id1=lcor(i,nelem)
   id2=lcor(i,maxelp)
   lcor(i,nelem)=id2
   lcor(i,maxelp)=id1
 10   CONTINUE
id1=idiv(nelem)
id2=idiv(maxelp)
idiv(nelem)=id2
idiv(maxelp)=id1
DO 20 k=1,9
   id1=id(k,nelem)
   id2=id(k,maxelp)
   id(k,nelem)=id2
   id(k,maxelp)=id1
 20   CONTINUE
RETURN
END
