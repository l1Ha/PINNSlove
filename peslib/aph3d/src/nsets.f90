SUBROUTINE nsets(n,ivalue,incx,iy,incy)
!
!   $RCSfile: nsets.f,v $ $Revision: 1.15 $
!   $Date: 89/11/16 11:03:15 $
!   $State: Stable $
!

!            P U R P O S E   O F    S U B R O U T I N E


!            I N P U T    A R G U M E N T S
!     n
!     ivalue
!     incx
!     iy
!     incy


IMPLICIT NONE



!            I N T E G E R S
INTEGER i, n, ivalue, incx, iy, incy


!            D I M E N S I O N S
DIMENSION iy(1)

!            C O M M O N S

DO 1 i = 1, n, incy
   iy(i) = ivalue
    1 CONTINUE

RETURN
END
