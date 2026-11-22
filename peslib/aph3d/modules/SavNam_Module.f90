MODULE SavNam_Module
   USE Numeric_Kinds_Module
   IMPLICIT NONE
   SAVE
   CHARACTER(LEN=70) ndir     !
   INTEGER(KIND=IW_Kind) :: IthSav=  100  !
   INTEGER(KIND=IW_Kind) LenDir  !
   NAMELIST/savnam/ ithsav, lendir, ndir
END MODULE SavNam_Module
