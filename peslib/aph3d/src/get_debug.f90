SUBROUTINE get_debug(procname_inp,debug)
USE FileNames_Module
USE FileUnits_Module
USE FileUnits_Common_Module
IMPLICIT NONE
CHARACTER(LEN=*) ::  procname_inp
LOGICAL       ::  debug

CHARACTER(LEN=1) , PARAMETER  :: per = ':'
CHARACTER(LEN=21), PARAMETER :: procname='get_debug'
CHARACTER(LEN=10), ALLOCATABLE:: subr(:)
CHARACTER(LEN=12)             :: InFormat
CHARACTER(LEN=10)             :: name1
CHARACTER(LEN=10)             :: name2
CHARACTER(LEN=30)             :: text
CHARACTER(LEN=30)             :: dummy
LOGICAL, ALLOCATABLE     :: debugv(:)
INTEGER i,nsubs,j
LOGICAL firstcall,isthere,metadebug

DATA firstcall /.true./
DATA metadebug /.true./
SAVE firstcall, debugv,subr

RETURN

IF(firstcall)THEN
    INQUIRE(file=OutDIR(1:LEN(TRIM(OutDIR)))//Debug_name,EXIST=isthere)
    IF(isthere)THEN
        OPEN(Unit= Dbug_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//Debug_name, status='old')
    ELSE
        CALL errormsg(Msg_Unit,procname,' Debug_name file does not exist.')
    ENDIF
    firstcall = .false.
    InFormat = '(A30,A40)'
    READ(Dbug_Unit,*) dummy                ! Comment
    READ(Dbug_Unit,*) text    ! Number of SUBROUTINEs
    WRITE(Msg_Unit,*) text
    READ(text,*) nsubs
    ALLOCATE (subr(nsubs))
    ALLOCATE (debugv(nsubs))
    InFormat = '(A10,A1,A40)'
    DO i=1,nsubs
        READ(Dbug_Unit,*) text
        !WRITE(Msg_Unit,*) text
        DO j = 1,len_trim(text)
            !WRITE(Msg_Unit,*) text(j:j)
            IF(text(j:j)==per)THEN
                subr(i) = text(1:j-1)
                READ(text(j+1:len_trim(text)),*) debugv(i)
                !WRITE(Msg_Unit,*) i,subr(i),'->',text,debugv(i)
                exit
            ENDIF
        ENDDO
        !READ(text,*) debugv(i)
    ENDDO
    CLOSE(Dbug_Unit)
ENDIF
!
! Locate the routine name
!
DO i=1, nsubs
    name1 = subr(i)
    name1 = TRIM(name1)   !name1(:Len_Trim(name1)) // char(0) !GregParker TEMPMOD
    name2 = TRIM(procname_inp)  !procname_inp(:Len_Trim(procname_inp)) // char(0) !GregParker TEMPMOD
    IF(name1 == name2)THEN
        debug = debugv(i)
        IF(metadebug)THEN
            IF(debug)THEN
                WRITE(Msg_Unit,'(A,A,A)') 'Debug for routine ',TRIM(procname_inp),' is ON'
            ELSE
                WRITE(Msg_Unit,'(A,A,A)') 'Debug for routine ',TRIM(procname_inp),' is OFF'
            ENDIF
        ENDIF
    ENDIF
ENDDO

RETURN
END
