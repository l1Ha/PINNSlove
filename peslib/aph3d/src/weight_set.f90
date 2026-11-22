SUBROUTINE weight_set(nang,nbasiss,philast,weight)
USE Numeric_Kinds_Module
USE FileUnits_Module
IMPLICIT NONE
INTEGER :: iang, ibasis, nang, nbasiss
REAL(Kind=WP_Kind) :: philast (nang, nbasiss), weight (nang)
DO 2 ibasis = 1, nbasiss
   DO 1 iang = 1, nang
      philast (iang, ibasis) = philast (iang, ibasis) * weight (iang)
    1    ENDDO
    2 ENDDO
RETURN
ENDSUBROUTINE weight_set
