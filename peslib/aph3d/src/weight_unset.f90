SUBROUTINE weight_unset(nang,nbasiss,philast,weight)
USE Numeric_Kinds_Module
USE FileUnits_Module
IMPLICIT NONE
INTEGER iang, ibasis, nang, nbasiss
REAL(Kind=WP_Kind) philast(nang,nbasiss),weight(nang)
DO ibasis=1,nbasiss
   DO iang=1,nang
      philast(iang,ibasis)=philast(iang,ibasis)/weight(iang)
   ENDDO
ENDDO
RETURN
ENDSUBROUTINE weight_unset
