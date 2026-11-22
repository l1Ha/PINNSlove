
!-------------------------------------------------------------------
! SUBROUTINE: warning
!
! package   : CID
!
! Language  : Fortran 90
!
! author    : F. Colavecchia (flavioc@lanl.gov)
!
! date      : Apr 02       version: 0.0
! revision  :              version:
!
! purpose   :  WRITEs a warning message.
!
! input     :    unit -> output file unit
!                procname -> name of the calling routine
!                message-> warning message
!
! output    :
!
! modules   :
!
! common    :
!
! notes     : usually procname is CHARACTER(LEN=10).
!
!-------------------------------------------------------------------
SUBROUTINE warning(unit,procname,message)
USE FileUnits_Module
IMPLICIT NONE
INTEGER unit
CHARACTER(LEN=*) procname
CHARACTER(LEN=*) message

WRITE(unit,*) "WARNING ",procname,": ",message
CALL flush(unit)
RETURN
END
