SUBROUTINE WriteVec(avec, n, outunit)
USE Numeric_Kinds_Module
USE FileUnits_Asymptotic_Module
IMPLICIT NONE
INTEGER n, outunit
REAL(Kind=WP_Kind) avec(n)
CHARACTER(LEN=21) end_label
CHARACTER(LEN=21) start_label
CHARACTER(LEN=21), PARAMETER:: ProcName='WriteVec'
CHARACTER(LEN=6) Print_Flag

CALL PoptAsy(ProcName, Print_Flag)
!
! Write the Vector.
!
start_label = 'Start Matrix         '
end_label =   'End Matrix           '
WRITE(outunit)start_label, avec, end_label

RETURN
ENDSUBROUTINE WriteVec
