FUNCTION delta (i, j)
USE Numeric_Kinds_Module
USE FileUnits_Module
IMPLICIT NONE
REAL(Kind=WP_Kind) :: delta
INTEGER :: i, j
IF(i==j)THEN
   delta = 1.d0
ELSE
   delta = 0.d0
ENDIF
RETURN
ENDFUNCTION delta
