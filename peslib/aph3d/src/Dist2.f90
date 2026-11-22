REAL(Kind=WP_Kind) FUNCTION Dist2(A,B,C)
USE Numeric_Kinds_Module
IMPLICIT NONE
REAL(Kind=WP_Kind) A, B, C
DIST2=B*B + (A + C)*(A + C)
RETURN
END
