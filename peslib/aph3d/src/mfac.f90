SUBROUTINE mfac (k1, k2, ffac)
!

!         P U R P O S E   O F    S U B R O U T I N E
! computes factorial(k1)/factorial(k2)
! assumes k1=n-l-1, k2=n+l+1

!         I N P U T    A R G U M E N T S
!  k1 = n - l - 1
!  k2 = n + l + 1

!         O U T P U T    A R G U M E N T S
! ffac = factorial(k1)/factorial(k2)
USE Numeric_Kinds_Module
USE FileUnits_Module
USE Numbers_Module
IMPLICIT NONE
REAL(Kind=WP_Kind)  ffac
INTEGER i, k1, k2

ffac = one
DO i = k1+1, k2
   ffac = ffac*i
ENDDO
ffac = one/ffac

RETURN
END
