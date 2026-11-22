REAL(Kind=WP_Kind) FUNCTION TimeInfo(Location, Calling)
USE Numeric_Kinds_Module
USE fileunits_Module
IMPLICIT NONE
!-------------------------------------------------------------------
! function  : TimeInfo
!
! package   : CID
!
! Language  : Fortran 90
!
! author    : G. A. Parker    (parker@nhn.ou.edu)
!             F. Colavecchia  (flavioc@lanl.gov)
!
!
! date      :              version: 0.0         G. A. Parker    (parker@nhn.ou.edu)
! revision  : Apr-2002     version: 0.1         F. Colavecchia  (flavioc@lanl.gov)
!
! purpose   :  Timing routine
!
! input     :    Location-> In which routine we are now
!                Calling -> Which routine is calling Location
!
! output    :    TimeInfo-> seconds spent since the first CALL
!
! modules   :
!             files -> file and units
!
! common    :   none
!
! functions : secnds
!
! notes     :  secnds is system-dependent
!
!
!-------------------------------------------------------------------

SAVE
LOGICAL,      PARAMETER:: DEBUG=.FALSE.
CHARACTER(LEN=1),  Parameter:: CTab ='#'
CHARACTER(LEN=15), PARAMETER:: Blank='               '
CHARACTER(LEN=*) Location
CHARACTER(LEN=30)  pLocation
CHARACTER(LEN=*) Calling
CHARACTER(LEN=8)  Date
CHARACTER(LEN=10) Time
CHARACTER(LEN=12) Now
CHARACTER(LEN=20) LastCall
REAL(KIND=WP_Kind) InitialTime, TimeNow
REAL(Kind=WP_Kind) :: DiffTime=0.0
REAL(Kind=SP_Kind) :: Secnds, TimeLast=0.0      ! These variables must be single precision
INTEGER ElapsedHours, ElapsedMinutes
EXTERNAL Secnds


CALL Date_And_Time(Date,Time)
Now=Time(1:2)//':'//Time(3:4)//':'//Time(5:6)//Time(7:10)

IF(DEBUG)THEN
    WRITE(Dbug_Unit,*)'Time: ',Now,' In ',Location,' Calling ',Calling
    WRITE(Out_Unit,*)'Time: ',Now,' In ',Location,' Calling ',Calling
ENDIF

IF(TimeLast==0.0)THEN
    TimeLast=TimeNow
    LastCall=Blank
ENDIF

!
! Do not print anything, just get time
!
IF(Calling=='')THEN
    TimeInfo = Secnds(TimeLast)
    RETURN
ENDIF

TimeNow=Secnds(TimeLast)
TimeInfo = TimeNow
DiffTime=TimeNow-DiffTime

IF(Location=='Main_Start')THEN
    TimeLast=TimeNow
    TimeNow=Secnds(TimeLast)
    InitialTime=TimeNow
    WRITE(Msg_Unit,*) '-------------------------------'
    WRITE(Msg_Unit,18)Date(5:6),Date(7:8),Date(1:4)
    WRITE(Msg_Unit,*)'Initial Time  =',Now
    WRITE(Msg_Unit,*)'Initial Second=',InitialTime
    WRITE(Msg_Unit,*) '-------------------------------'
    WRITE(Out_Unit,18)Date(5:6),Date(7:8),Date(1:4)
    WRITE(Out_Unit,*)'Initial Time  =',Now
    WRITE(Out_Unit,*)'Initial Second=',InitialTime
    WRITE(Out_Unit,*) '-------------------------------'
    WRITE(Time_Unit,18)Date(5:6),Date(7:8),Date(1:4)
    WRITE(Time_Unit,*)'Initial Time  =',Now
    WRITE(Time_Unit,*)'Initial Second=',InitialTime
    WRITE(Time_Unit,*)'  Time (sec)    Routine     Time in routine    Calling'
    WRITE(Time_Unit,*)'--------------------------------------------------------'
ENDIF
plocation = location
WRITE(Time_Unit,'(es15.7,4x,a1,4x,A,4x,a1,a1,4x,A)')TimeNow, CTab, TRIM(pLocation), CTab, CTab, Calling
IF(len(Calling)>=5)THEN
   IF(Calling(1:5)=='  End')THEN
       WRITE(Time_Unit,'(4x,A,es15.7,A,es15.7,A,A)')"TimeNow=", TimeNow, "  DiffTime=", DiffTime, "  Time in ", TRIM(pLocation)
       WRITE(Time_Unit,*)
   ENDIF
ENDIF

IF(Location=='Main_Stop')THEN
    ElapsedHours   = (TimeNow-InitialTime)/3600
    ElapsedMinutes = ((TimeNow-InitialTime)-ElapsedHours*3600)/60
    WRITE(Msg_Unit,*)'Final Second      =',TimeNow
    WRITE(Msg_Unit,*)'Final Time        =',Now
    WRITE(Msg_Unit,*)'Total Elapsed Time=',TimeNow-InitialTime
    WRITE(Msg_Unit,'(1x,a19,i5,a4,i5,a7,i7,a3)')'Total Walltime    =',ElapsedHours,' hs ',ElapsedMinutes,' min = ', &
                                        int(TimeNow-InitialTime),' s.'
    WRITE(Msg_Unit,18)Date(5:6),Date(7:8),Date(1:4)
    WRITE(Msg_Unit,19)
    WRITE(Out_Unit,*)'Final Second      =',TimeNow
    WRITE(Out_Unit,*)'Final Time        =',Now
    WRITE(Out_Unit,*)'Total Elapsed Time=',TimeNow-InitialTime
    WRITE(Out_Unit,'(1x,a19,i5,a4,i5,a7,i7,a3)')'Total Walltime    =',ElapsedHours,' hs ',ElapsedMinutes,' min = ', &
                                        int(TimeNow-InitialTime),' s.'
    WRITE(Out_Unit,18)Date(5:6),Date(7:8),Date(1:4)
    WRITE(Out_Unit,19)
    WRITE(Time_Unit,*)'Final Second      =',TimeNow
    WRITE(Time_Unit,*)'Final Time        =',Now
    WRITE(Time_Unit,*)'Total Elapsed Time=',TimeNow-InitialTime
    WRITE(Time_Unit,'(1x,a19,i5,a4,i5,a7,i7,a3)')'Total Walltime    =',ElapsedHours,' hs ',ElapsedMinutes,' min = ', &
                                        int(TimeNow-InitialTime),' s.'
    WRITE(Time_Unit,18)Date(5:6),Date(7:8),Date(1:4)
    WRITE(Time_Unit,19)
ENDIF


DiffTime=TimeNow
LastCall=Blank
LastCall=Calling

CALL Flush(Time_Unit)

RETURN
17 FORMAT(1x,'Starting')
18 FORMAT(1x,'Date: ',a2,'/',a2,'/',a4)
19 FORMAT(1x,'Finished')
END FUNCTION TimeInfo
