SUBROUTINE NameList_Default(Out_Unit, NameList_Name)
IMPLICIT NONE
CHARACTER(*) NameList_Name
INTEGER Out_Unit
LOGICAL, PARAMETER :: PRT_NameList=.FALSE.
    
IF(PRT_NameList)THEN
   WRITE(Out_Unit,*)'Namelist '//TRIM(NameList_Name)//' using default values'
ENDIF

RETURN
END SUBROUTINE NameList_Default