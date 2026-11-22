subroutine WRITE_real(unit,comment,dat)
USE Numeric_Kinds_Module
IMPLICIT NONE
integer unit
CHARACTER(LEN=*) comment
CHARACTER(LEN=10) text,text2
REAL(Kind=WP_Kind) dat

WRITE(text2,'(1f18.8)') dat
text = comment(:len_trim(comment))//text2
WRITE(unit,'(a)') text(:len_trim(text))
call flush(unit)

return
end
