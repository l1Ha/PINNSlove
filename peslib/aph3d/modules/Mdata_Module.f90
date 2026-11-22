MODULE Mdata_Module
   USE Numeric_Kinds_Module
   IMPLICIT NONE
   SAVE
   REAL(Kind=WP_Kind), PARAMETER:: rtol=1.d-7  !
   INTEGER(KIND=IW_Kind) :: nlt     = 17                     !
   INTEGER(KIND=IW_Kind) :: nlc(3)  = 8                 !
   INTEGER(KIND=IW_Kind) :: nq                       !
   INTEGER(KIND=IW_Kind) :: nmode                    !
   INTEGER(KIND=IW_Kind) :: nodmax  = 19000                 !
   INTEGER(KIND=IW_Kind) :: nblock                   !
   INTEGER(KIND=IW_Kind) :: maxj                     !
   INTEGER(KIND=IW_Kind) :: nlancz  = 2                 !
END MODULE Mdata_Module
