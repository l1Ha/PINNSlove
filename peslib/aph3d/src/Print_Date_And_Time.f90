SUBROUTINE Print_Date_And_Time(Today, Hour, CurZone, DtValues, Out_Unit)
IMPLICIT NONE
CHARACTER(LEN=10) hour
CHARACTER(LEN=8) today
CHARACTER(LEN=5) curzone
CHARACTER(LEN=21), PARAMETER:: ProcName='Print_Date_And_Time'
CHARACTER(LEN=6) Print_Flag
INTEGER Out_Unit, dtvalues(8)
LOGICAL PRT_Timing/.FALSE./

CALL PoptAsy(ProcName, Print_Flag)
IF(Print_Flag=='Full')PRT_Timing=.True.
IF(PRT_Timing)WRITE(Out_Unit,1)today, hour, curzone, dtvalues

RETURN
1  FORMAT(1x,'Date: ',a10,/, 1x,'Time: ',a10,/, 1x,'Zone: ',a10,/, &
          1x,'year: ',i4,/, 1x,'Month:',i3,/, 1x,'Day:',i3,/,    &
          1x,'Time Difference:',i6,/, &
          1x,'Hour:',i3,/, 1x,'Minute:',i3,/, 1x,'Second:',i3,/, 1x,'Millisecond:',i4,/,/)
END SUBROUTINE Print_Date_And_Time