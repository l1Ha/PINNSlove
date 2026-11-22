SUBROUTINE corner (iel, lcor, nnc, jmax, iskp, jskp, n, idiv, ibox, jbox, dthe, dchi, dang, mxelmnt, icor, jcor)
USE Numeric_Kinds_Module
USE FileUnits_Module
!
!     $RCSfile: corner.f,v $ $Revision: 1.14 $
!     $Date: 89/10/18 14:17:58 $
!     $State: Stable $
!


!     I N P U T    A R G U M E N T S
!
!     iel
!     lcor
!     nnc
!     jmax
!     iskp
!     jskp
!     n
!     idiv
!     ibox
!     jbox
!     dchi
!     dang
!     mxelmnt
!     icor
!     jcor


IMPLICIT NONE

!     I N T E G E R S

INTEGER ibox, icor, idiv, iel, iskp, jbox, jcor, jmax, jskp, lcor, mxelmnt, n, nnc

!     R E A L S
REAL(Kind=WP_Kind)  dang, dchi, dthe

!     D I M E N S I O N S

DIMENSION lcor(2,mxelmnt), nnc(3), jskp(3), idiv(mxelmnt), dchi(3)

!     I N T R I N S I C    F U N C T I O N S
INTRINSIC sqrt

!-----------------------------------------------------------------------
!     find the corner of the current element; which region of chi
!     it is in, the number of nodes between its turned on nodes, and its
!     corner to corner distance
!-----------------------------------------------------------------------
icor=lcor(1,iel)
jcor=lcor(2,iel)
IF(jcor<nnc(1))THEN
   n=1
ELSEIF (jcor<(nnc(1)+nnc(2)-1))THEN
   n=2
ELSEIF (jcor<jmax)THEN
   n=3
ELSE
   WRITE(Out_unit,*) '*****error*****: jcor>=jmax !','jcor=',jcor,  ' jmax= ',jmax
   STOP 'corner'
ENDIF
!-----------------------------------------------------------------------
! ibox and jbox are the number of possible nodes between two turned on n
! plus one for one of the turned on nodes
!-----------------------------------------------------------------------
ibox=iskp/(2**(idiv(iel)+1))
jbox=jskp(n)/(2**(idiv(iel)+1))
dang=sqrt((ibox*dthe)**2+(dchi(n)*jbox)**2)
RETURN
END
