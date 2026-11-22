
    subroutine WRITE_int(unit,comment,dat)
    IMPLICIT NONE
    integer unit,dat
    CHARACTER(LEN=*) comment
    CHARACTER(LEN=10) text,text2

    WRITE(text2,'(i5)') dat
    text = comment(:len_trim(comment))//text2
    WRITE(unit,'(a)') text(:len_trim(text))
    call flush(unit)

    return
    end
