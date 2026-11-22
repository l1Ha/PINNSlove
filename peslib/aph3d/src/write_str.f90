
    subroutine WRITE_str(unit,comment,text)
    IMPLICIT NONE
    integer unit
    CHARACTER(LEN=*) comment
    CHARACTER(LEN=*) text
    WRITE(unit,'(a)') comment(:len_trim(comment))//text(:len_trim(text))
    call flush(unit)
    return
    end
