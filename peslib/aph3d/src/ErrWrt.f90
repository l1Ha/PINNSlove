SUBROUTINE ErrWrt (ErrMsg)
USE FileUnits_Module
!
!-------------------------------- ErrWrt -------------------------------
!
! >>> ErrWrt   -- WRITEs a message to the standard error output file
!
! *** Calling Sequence:   CALL ErrWrt ( ErrMsg )
!
! === Latest Revision:    25 May 1988  (typed in)
!
! --- on entry the calling routine supplies --
!
!     ErrMsg   -- CHARACTER string to be written to standard error
!
! --- on exit the routine returns --
!
!     nothing
!
!----------------------------- local declarations ----------------------
!
!
CHARACTER(LEN=*) ::    ErrMsg
!
! * * * * *                !  INTEGER globals
CHARACTER(LEN=8  ) :: IGlobal
! * * * * *              !    output line
CHARACTER(LEN=80    ) :: Line
! * * * * *               !   output FORMAT
CHARACTER(LEN=8  ) :: MyForm
!
! * * * * *            !   standard error unit number
INTEGER ::        StdErr
! * * * * *            !   right index function
INTEGER ::        IndexR
! * * * * *           !    # CHARACTERs to output
INTEGER ::        NCOut
!
!---------------------------- execution begins here --------------------
!
IF(IGlobal ('get', 'StdErr', StdErr) /='OK') RETURN
!
Line   = ErrMsg
NCOut     = max (1, IndexR (Line) )
!
MyForm    = '(a!!)'
WRITE(MyForm (3:4) , '(i2)') NCOut
!
WRITE(StdErr, MyForm) Line (1:NCOut)
!
RETURN
ENDSUBROUTINE ErrWrt
