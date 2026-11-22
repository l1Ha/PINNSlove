SUBROUTINE SUBV(X, Y, C, N)
USE Numeric_Kinds_Module
!
!   $RCSfile: subv.f,v $ $Revision: 1.3 $
!   $Date: 89/07/28 09:57:46 $
!   $State: Stable $
!
!
! **********************************************************************
!
!     PURPOSE - (VER = 1 OR 2)
!
!         COMPUTES X=X-C*Y.
!
!
!     INPUT PARAMETERS -
!
!         X AND Y = N-VECTORS PARTICIPATING IN THE COMPUTATION.
!         C       = REAL CONSTANT.
!         N       = DIMENSION OF THE VECTORS (NOT NECESSARILY OF THE
!                   PROBLEM).
!
!
!     OUTPUT PARAMETERS -
!
!         X = UPDATED X.
!
! **********************************************************************
!
REAL(Kind=WP_Kind) C, X(1), Y(1)
INTEGER   I, N
!
DO 10 I = 1, N
  X(I) = X(I) - C * Y(I)
10      CONTINUE
!
RETURN
END
