SUBROUTINE modcor (xx)
USE Numeric_Kinds_Module
!
!   $RCSfile: modcor.f,v $ $Revision: 1.15 $
!   $Date: 89/11/16 11:03:14 $
!   $State: Stable $
!

!            P U R P O S E   O F    S U B R O U T I N E


!            I N P U T    A R G U M E N T S
!    xx


IMPLICIT NONE



!            I N T E G E R S
INTEGER i, n2

!            R E A L S
REAL(Kind=WP_Kind) xx, fac

!            D I M E N S I O N S
DIMENSION xx(2,9), fac(8)


!     node 9 for serendipity elements
!
DO 70 i=1,4
fac(i)=.25d0
fac(i+4)=0.0d0
   70 CONTINUE
DO 80 i=1,4
n2=i+1
IF(n2>4) n2=1
fac(i)=fac(i)-.25d0
fac(n2)=fac(n2)-.25d0
fac(i+4)=.5d0
   80 CONTINUE
xx(1,9)=0.0d0
xx(2,9)=0.0d0
DO 90 i=1,8
xx(1,9)=xx(1,9)+xx(1,i)*fac(i)
xx(2,9)=xx(2,9)+xx(2,i)*fac(i)
   90 CONTINUE
RETURN
END
