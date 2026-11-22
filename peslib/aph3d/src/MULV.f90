SUBROUTINE MULV(X, Y, C, N)
USE Numeric_Kinds_Module
!
!   $RCSfile: mulv.f,v $ $Revision: 1.3 $
!   $Date: 89/07/28 09:57:24 $
!   $State: Stable $
!
!
! **********************************************************************
!
!     PURPOSE - (VER = 1 OR 2)
!
!         COMPUTES X=C*Y.
!
!
!     INPUT PARAMETERS -
!
!         Y = INCOMING VECTOR.
!         C = REAL CONSTANT.
!         N = DIMENSION OF THE VECTOR (NOT NECESSARILY OF THE PROBLEM).
!
!
!     OUTPUT PARAMETERS -
!
!         X = RESULTING VECTOR, C*Y.
!
! **********************************************************************
!
REAL(Kind=WP_Kind)  C, X(1), Y(1)
INTEGER   I, N
!
DO 10 I = 1, N
  X(I) = C * Y(I)
10      CONTINUE
!
RETURN
END
