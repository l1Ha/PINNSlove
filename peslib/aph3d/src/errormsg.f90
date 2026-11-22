SUBROUTINE ErrorMsg(unit,procname,message)
USE FileNames_Module
IMPLICIT NONE
INTEGER unit
CHARACTER(LEN=*) procname
CHARACTER(LEN=*) message
 
WRITE(unit,*)'ERROR in ', procname, ": ", message
CALL flush(Out_Unit)
STOP 'ErrorMsg'
ENDSUBROUTINE ErrorMsg
