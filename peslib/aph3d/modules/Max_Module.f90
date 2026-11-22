MODULE Max_Module
   USE Numeric_Kinds_Module
   IMPLICIT NONE
   SAVE
   INTEGER(KIND=IW_Kind), PARAMETER:: maxang=336                  !
   INTEGER(KIND=IW_Kind), PARAMETER:: maxangv=maxang*(maxang+1)/2 !
END MODULE Max_Module
