SUBROUTINE aph_to_delves
!
!
!            P U R P O S E   O F    S U B R O U T I N E
! This is the main program for the aph to delves transformation.
!
!            I N P U T    A R G U M E N T S
!  NONE
!
!            O U T P U T    A R G U M E N T S
!  NONE
!
!this routine calls:
!     date_and_time,title,constant,readaph,openaph,aph3d
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

USE FileUnits_Module
USE Numeric_Kinds_Module
USE Narran_Module
USE convrsns_Module
!USE Gaussb_Module
USE GaussQuady_Module
USE NState_Module
USE Parms_Module
IMPLICIT NONE
CHARACTER(LEN=10) hour
CHARACTER(LEN=8) today
CHARACTER(LEN=5) curzone
INTEGER dtvalues(8)

REAL(Kind=WP_Kind) w(nstate,nstate)
REAL(Kind=WP_Kind) vecnew(nstate,nstate)
REAL(Kind=WP_Kind) tstore(nstate,nstate)
REAL(Kind=WP_Kind) vecold(nstate,nstate)
REAL(Kind=WP_Kind) xsq(nstate)
REAL(Kind=WP_Kind) xk(nstate)

EXTERNAL title, constant, readaph, aph3d,  openaph

!----------------------------------------------------------------------
! Get primary input data.
!----------------------------------------------------------------------

WRITE(Out_Unit,*)'Called    APH_To_Delves'
CLOSE(Out_Unit)
OPEN(Unit=Out_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'Output/APH_To_Delves.txt',Form='FORMATTED')

!-----------------------------------------------------------------------
! Determine the time and todays date.
!-----------------------------------------------------------------------
CALL date_and_time(today, hour, curzone, dtvalues)
CALL Print_Date_And_Time(Today, Hour, CurZone, DtValues, Out_Unit)
!-----------------------------------------------------------------------
!   WRITE out the title page.
!-----------------------------------------------------------------------
CALL title

!-----------------------------------------------------------------------
! Establish all conversion factors and physical constants
!-----------------------------------------------------------------------
CALL constant

!-----------------------------------------------------------------------
! Read in all input data
!-----------------------------------------------------------------------
CALL readaph
IF(.not.Allocated(xpth))ALLOCATE(xpth(maxhermt,narran), wpth(maxhermt,narran))
!-----------------------------------------------------------------------
! Open files for input and output of data.
!-----------------------------------------------------------------------
CALL openaph

!-----------------------------------------------------------------------
! Propagate the set of coupled radial equations to determine transition
! probabilities and S-matrices.
!-----------------------------------------------------------------------
REWIND elm_unit
REWIND nmode_unit
CALL aph3d (w, vecnew, tstore, vecold, xsq, xk, nstate)
!
!DEALLOCATE(xpth,wpth)
CALL date_and_time(today, hour, curzone, dtvalues)
!-------------------------------------------------------------
!   WRITE the time and date to Out_Unit.
!-------------------------------------------------------------
CALL Print_Date_And_Time(Today, Hour, CurZone, DtValues, Out_Unit)

WRITE(Out_Unit,*)'Completed APH_To_Delves'
WRITE(Out_Unit,*)

CLOSE(Out_Unit)

OPEN(Unit=Out_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'Output/BasisForScatt.txt',Form='FORMATTED',ACCESS='APPEND')
WRITE(Out_Unit,*)'Completed APH_To_Delves'
WRITE(Out_Unit,*)

RETURN
ENDSUBROUTINE APH_to_Delves
