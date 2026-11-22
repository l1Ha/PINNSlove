SUBROUTINE PoptAsy(Subr, Print_Flag)
USE POptAsy_Module
USE FileUnits_Asymptotic_Module
IMPLICIT NONE
! Author: Gregory A. Parker, Department of Physics and Astronomy, University of Oklahoma
!
!  This routine sets the following parameters for printing options:
!  None           No matrix output
!  Little         The first row of the each matrix
!  Medium         A small block of the each matrix
!  Full           All possible output and the entire matrix
!
!  Required Input <=====
!     Subr        ProcName as given in the routine
!
!  On return  =====>
!     Print_Flag  Desire print level (None, Little, Medium, Full)
!
INTEGER Index
CHARACTER(LEN=*) Subr, Print_Flag

return !tmpmod gregparker

IF(AllFull)THEN
   Print_Flag='Full'
   RETURN
ELSEIF(AllMedium)THEN
   Print_Flag='Medium'
   RETURN
ELSEIF(AllLittle)THEN
   Print_Flag='Little'
   RETURN
ELSEIF(ALLNone)THEN
   Print_Flag='None'
   RETURN
ENDIF

DO index=1,NSubs
   IF(TRIM(Subr)==TRIM(Subs(Index)))THEN
      Print_Flag=Print_Level(Index)
      RETURN
   ENDIF
ENDDO

WRITE(Msg_Unit,*)"Error in POptAsy:", Subr,Print_Flag
WRITE(Out_Unit,*)"Error in POptAsy:", Subr,Print_Flag
STOP 'POptAsy'
ENDSUBROUTINE POptAsy
