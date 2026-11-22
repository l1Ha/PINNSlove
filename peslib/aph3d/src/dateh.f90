SUBROUTINE Dateh(Today)
USE Numeric_Kinds_Module
USE FileUnits_Module
IMPLICIT NONE
INTEGER           :: DtValues(8)
CHARACTER(LEN= 8) :: Today
CHARACTER(LEN=10) :: Hour
CHARACTER(LEN= 5) :: CurZone

CALL Date_And_Time(DATE=Today, TIME=Hour, ZONE=CurZone, VALUES=DtValues)
CALL Print_Date_And_Time(Today, Hour, CurZone, DtValues, Out_Unit)

RETURN
END
