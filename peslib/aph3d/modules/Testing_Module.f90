MODULE Testing_Module
   USE Numeric_Kinds_Module
   IMPLICIT NONE
   SAVE
   CHARACTER(LEN=3) :: test='nul'     !
   NAMELIST/Testing/ test
END MODULE Testing_Module
