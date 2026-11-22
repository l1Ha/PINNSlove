SUBROUTINE repnames (irrep, group, repname)
USE FileUnits_Module
IMPLICIT NONE
INTEGER :: irrep

CHARACTER(LEN=3) :: group, repname
IF(group=='C1 ')THEN
   IF(irrep==1)THEN
      repname = 'Axx'
   ELSE
      WRITE(Out_Unit, * ) 'irrep not correct: irrep=', irrep
   ENDIF
ELSEIF (group=='C2 ')THEN
   IF(irrep==1)THEN
      repname = 'Axx'
   ELSEIF (irrep==2)THEN
      repname = 'Bxx'
   ELSE
      WRITE(Out_Unit, * ) 'irrep not correct: irrep=', irrep
   ENDIF
ELSEIF (group=='C1v')THEN
   IF(irrep==1)THEN
      repname = 'A1x'
   ELSEIF (irrep==2)THEN
      repname = 'A2x'
   ELSE
      WRITE(Out_Unit, * ) 'irrep not correct: irrep=', irrep
   ENDIF
ELSEIF (group=='C2v')THEN
   IF(irrep==1)THEN
      repname = 'A1x'
   ELSEIF (irrep==2)THEN
      repname = 'A2x'
   ELSEIF (irrep==3)THEN
      repname = 'B1x'
   ELSEIF (irrep==4)THEN
      repname = 'B2x'
   ELSE
      WRITE(Out_Unit, * ) 'irrep not correct: irrep=', irrep
   ENDIF
ELSEIF (group=='C3 ')THEN
   IF(irrep==1)THEN
      repname = 'Axx'
   ELSEIF (irrep==2)THEN
      repname = 'Exa'
   ELSEIF (irrep==3)THEN
      repname = 'Exb'
   ELSE
      WRITE(Out_Unit, * ) 'irrep not correct: irrep=', irrep
   ENDIF
ELSEIF (group=='C3v')THEN
   IF(irrep==1)THEN
      repname = 'A1x'
   ELSEIF (irrep==2)THEN
      repname = 'A2x'
   ELSEIF (irrep==3)THEN
      repname = 'Exa'
   ELSEIF (irrep==4)THEN
      repname = 'Exb'
   ELSE
      WRITE(Out_Unit, * ) 'irrep not correct: irrep=', irrep
   ENDIF
ELSEIF (group=='C6 ')THEN
   IF(irrep==1)THEN
      repname = 'Axx'
   ELSEIF (irrep==2)THEN
      repname = 'Bxx'
   ELSEIF (irrep==3)THEN
      repname = 'E1a'
   ELSEIF (irrep==4)THEN
      repname = 'E1b'
   ELSEIF (irrep==5)THEN
      repname = 'E2a'
   ELSEIF (irrep==6)THEN
      repname = 'E2b'
   ELSE
      WRITE(Out_Unit, * ) 'irrep not correct: irrep=', irrep
   ENDIF
ELSEIF (group=='C6v')THEN
   IF(irrep==1)THEN
      repname = 'A1x'
   ELSEIF (irrep==2)THEN
      repname = 'A2x'
   ELSEIF (irrep==3)THEN
      repname = 'B1x'
   ELSEIF (irrep==4)THEN
      repname = 'B2x'
   ELSEIF (irrep==5)THEN
      repname = 'E1a'
   ELSEIF (irrep==6)THEN
      repname = 'E1b'
   ELSEIF (irrep==7)THEN
      repname = 'E2a'
   ELSEIF (irrep==8)THEN
      repname = 'E2b'
   ELSE
      WRITE(Out_Unit, * ) 'irrep not correct: irrep=', irrep
   ENDIF
ELSE
   WRITE(Out_Unit, * ) 'Incorret or no group: group=', group
   STOP 'symop'
ENDIF
RETURN
ENDSUBROUTINE repnames
