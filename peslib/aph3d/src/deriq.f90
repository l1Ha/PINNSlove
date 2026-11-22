SUBROUTINE deriq (nel, xx, b, det, r, s, x1bar, theta, phi)
USE FileUnits_Module
USE todim_Module
!
!   $RCSfile: deriq.f,v $ $Revision: 1.15 $
!   $Date: 89/11/16 11:02:57 $
!   $State: Stable $
!

!            P U R P O S E   O F    S U B R O U T I N E


!            I N P U T    A R G U M E N T S

!     nel
!     xx
!     b
!     det
!     r
!     s
!     x1bar
!     theta
!     phi


IMPLICIT NONE

!            I N T E G E R S
INTEGER nel, iintp, k,i, k2

!            R E A L S
REAL(Kind=WP_Kind) zero,one, dum
REAL(Kind=WP_Kind) xx, b, det, r, s, x1bar, theta, phi, h, p, xj, xji

!            D I M E N S I O N S
DIMENSION xx(2,1), b(4,1), h(9), p(2,9), xj(2,2), xji(2,2)

!            E X T E R N A L S
PARAMETER (zero=0.d0,one=1.d0)

!-----------------------------------------------------------------------
!     find interpolation functions and jacobian
!-----------------------------------------------------------------------
iintp=0
CALL funct2 (r, s, h, p, xj, det, xx, iintp)
!-----------------------------------------------------------------------
!     compute inverse of the jacobian matrix
!-----------------------------------------------------------------------
dum=one/det
xji(1,1)=xj(2,2)*dum
xji(1,2)=-xj(1,2)*dum
xji(2,1)=-xj(2,1)*dum
xji(2,2)=xj(1,1)*dum
!-----------------------------------------------------------------------
!     evaluate global derivative operator ( b-matrix )
!-----------------------------------------------------------------------
DO k=1,9
k2=k*2
b(1,k2-1)=zero
b(1,k2)=zero
b(2,k2-1)=zero
b(2,k2)=zero
DO i=1,2
b(1,k2-1)=b(1,k2-1)+xji(1,i)*p(i,k)
b(2,k2)=b(2,k2)+xji(2,i)*p(i,k)
ENDDO
b(3,k2)=b(1,k2-1)
b(3,k2-1)=b(2,k2)
ENDDO
x1bar=one
theta=zero
phi=zero
DO k=1,9
theta=theta+h(k)*xx(1,k)
phi=phi+h(k)*xx(2,k)
ENDDO
IF(x1bar>.00000001) GOTO 50
!-----------------------------------------------------------------------
!     for the case of zero radius equate radial to hoop strain
!-----------------------------------------------------------------------
DO k=1,18
   b(4,k)=b(1,k)
ENDDO
RETURN
!-----------------------------------------------------------------------
!     non-zero radius
!-----------------------------------------------------------------------
   50 dum=one
DO k=1,9
k2=k*2
b(4,k2)=zero
b(4,k2-1)=h(k)
ENDDO
!
RETURN
END
