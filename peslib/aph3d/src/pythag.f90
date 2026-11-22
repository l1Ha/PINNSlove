REAL(Kind=WP_Kind) FUNCTION PYTHAG(A,B)
USE Numeric_Kinds_Module
USE FileUnits_Module
REAL(Kind=WP_Kind) A, B, P,Q,R,S,T
!   BEGIN PROLOGUE  PYTHAG
!   REFER TO  EISDOC
!
!     Finds sqrt(A**2+B**2) without overflow or destructive underflow
!   ROUTINES CALLED  (NONE)
!   END PROLOGUE  PYTHAG
!
!     P = AMAX1(ABS(A),ABS(B))
P = MAX(ABS(A),ABS(B))
!     Q = AMIN1(ABS(A),ABS(B))
Q = MIN(ABS(A),ABS(B))
IF(Q == 0.0D0) GO TO 20
   10 CONTINUE
   R = (Q/P)**2
   T = 4.0D0 + R
   IF(T == 4.0D0) GO TO 20
   S = R/T
   P = P + 2.0D0*P*S
   Q = Q*S
GO TO 10
   20 PYTHAG = P
RETURN
END FUNCTION pythag
