SUBROUTINE center (xx)
USE Numeric_Kinds_Module
!
!   $RCSfile: center.f,v $ $Revision: 1.15 $
!   $Date: 89/11/16 11:02:51 $
!   $State: Stable $
!

!            P U R P O S E   O F    S U B R O U T I N E

! Make the coordinate of the sides nodes at the midpoint of the
! corner nodes.


!            I N P U T    A R G U M E N T S
!  xx


IMPLICIT NONE


!            I N T E G E R S
INTEGER i, n2

!            R E A L S
REAL(Kind=WP_Kind) xx

!            D I M E N S I O N S
DIMENSION xx(2, 1)

DO i=1,4
   n2=i+1
   IF(n2==5) n2=1
   xx(1,i+4)=.5d0*(xx(1,i)+xx(1,n2))
   xx(2,i+4)=.5d0*(xx(2,i)+xx(2,n2))
ENDDO
RETURN
END
