SUBROUTINE check_label(start_label, end_label, label)
USE FileUnits_Module
IMPLICIT NONE
!#include <Inc.unitno>
CHARACTER(LEN=11) label
CHARACTER(LEN=21) end_label
CHARACTER(LEN=21) start_label

!
! Check start label.
!
IF(start_label/='Start '//label)THEN
    WRITE(Msg_Unit,*)'ERROR: in routine check_label'
    WRITE(Msg_Unit,*)'Label = ',label
    WRITE(Msg_Unit,*)'Start_label = ', start_label
    WRITE(Msg_Unit,*)'End_label = ', end_label
    STOP 'check_label'
ENDIF

!
! Check end label.
!
IF(end_label/='End '//label)THEN
    WRITE(Msg_Unit,*)'ERROR: in routine check_label'
    WRITE(Msg_Unit,*)'Label = ',label
    WRITE(Msg_Unit,*)'Start_label = ', start_label
    WRITE(Msg_Unit,*)'End_label = ', end_label
    STOP 'check_label'
ENDIF

RETURN
ENDSUBROUTINE Check_Label
