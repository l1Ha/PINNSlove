SUBROUTINE doubpi(nelem, th, ch, maxelp, mxelmnt, mapr)
USE FileUnits_Module
USE Numbers_Module
!-----------------------------------------------------------------------
!     This routine doubles the mesh for symmetry=false; any point at
!     (theta,chi) is duplicated at (theta,chi+pi).
!-----------------------------------------------------------------------

!     I N P U T    A R G U M E N T S

!     nelem
!     th
!     ch
!     maxelp
!     mxelmnt
!     mapr


IMPLICIT NONE
INTEGER, INTENT(IN)               :: mxelmnt
INTEGER, INTENT(INOUT)            :: mapr(mxelmnt)
INTEGER, INTENT(INOUT)            :: nelem, maxelp
INTEGER                           :: iel, ielp, k
REAL(Kind=WP_Kind), INTENT(INOUT) :: th(9,mxelmnt), ch(9,mxelmnt)

IF(2*nelem>mxelmnt)THEN
   WRITE(Out_unit,*)'***error***: 2*nelem>mxelmnt; nelem=',nelem,  ' mxelmnt=',mxelmnt
   STOP 'doubpi'
ENDIF

ielp=nelem
DO iel=1,nelem
   ielp=ielp+1
   mapr(ielp)=iel
   DO k=1,9
      th(k,ielp)=th(k,iel)
      ch(k,ielp)=ch(k,iel)+pi
   ENDDO
ENDDO
maxelp=ielp
nelem=2*nelem
!
RETURN
END
