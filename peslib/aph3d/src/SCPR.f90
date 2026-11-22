REAL(Kind=WP_Kind) FUNCTION SCPR(X, Y, N)
USE Numeric_Kinds_Module
!
!   $RCSfile: scpr.f,v $ $Revision: 1.3 $
!   $Date: 89/07/28 09:57:39 $
!   $State: Stable $
!
!
! **********************************************************************
!
!     PURPOSE - (VER = 1 OR 2)
!
!         COMPUTES SCPR=X(T)*Y.
!
!
!     INPUT PARAMETERS -
!
!         X AND Y = VECTORS PARTICIPATING IN THE COMPUTATION.
!         N       = THE DIMENSION OF THE VECTORS (NOT NECESSARILY OF THE
!                   PROBLEM).
!
!
!     OUTPUT PARAMETERS -
!
!         SCPR = THE SCALARPRODUCT BETWEEN THE VECTORS.
!
! **********************************************************************
!
REAL(Kind=WP_Kind) X(1), Y(1)
INTEGER   I, N
!
SCPR = 0.0D0
DO 10 I = 1, N
  SCPR = SCPR + X(I) * Y(I)
10      CONTINUE
!
RETURN
END
