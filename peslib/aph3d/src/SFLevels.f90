SUBROUTINE SFLevels(In_File)
USE FileUnits_Module
USE Boundary_Module
IMPLICIT NONE
INTEGER dtvalues(8)
CHARACTER(LEN=21), PARAMETER:: ProcName='SFLevels'
CHARACTER(LEN=5 ) readx, ready
CHARACTER(LEN=8 ) cont
CHARACTER(LEN=10) today, hour, curzone
CHARACTER(LEN=*) In_File
LOGICAL there
INTEGER ieng, neng, ienergy, nstate, iunit, ounit, istart, nchanl, mxeng
INTEGER, PARAMETER:: mxstate=300
REAL, ALLOCATABLE::  engvals(:), evalues(:)
EXTERNAL wgdat

WRITE(Out_Unit,*)'Called    SFLevels'
CLOSE(Out_Unit)
OPEN(Unit=Out_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'Output/SFLevels.txt',Form='FORMATTED')
mxeng=Max_Sector+5
ALLOCATE(engvals(mxeng), evalues(mxeng*mxstate))

WRITE(Out_Unit,*)
WRITE(Out_Unit,*)
WRITE(Out_Unit,*)
WRITE(Out_Unit,*)'Entering:', ProcName

CALL Date_And_Time(today, hour, curzone, dtvalues)
CALL Print_Date_And_Time(Today, Hour, CurZone, DtValues, Out_Unit)

ieng = 1
istart = 1
iunit=21
ounit=91

INQUIRE(file=OutDIR(1:LEN(TRIM(OutDIR)))//In_File, exist=there)
IF(there)THEN
   OPEN(Unit=iunit,File=OutDIR(1:LEN(TRIM(OutDIR)))//In_File,status='old')
ELSE
   WRITE(Out_Unit,*)OutDIR(1:LEN(TRIM(OutDIR)))//In_File, 'No sfelevl file exists: Stopped in SFLevels'
   WRITE(Msg_Unit,*)OutDIR(1:LEN(TRIM(OutDIR)))//In_File, 'No sfelevl file exists: Stopped in SFLevels'
   STOP 'sfelevl'
ENDIF

30 READ(iunit,8,END=2)readx
IF(readx/='readx')THEN
   WRITE(Out_Unit,*)'readx is not correct'
   WRITE(Out_Unit,*)'Stopping in SFLevels'
   STOP 'sfelevl'
ENDIF

READ(iunit,9,END=2)neng
IF(neng/=1)THEN
   WRITE(Out_Unit,*)'neng/=1; neng=',neng
   STOP 'sfelevl'
ENDIF

READ(iunit,3,END=2)engvals(ieng)
READ(iunit,8,END=2)ready

IF(ready/='ready')THEN
   WRITE(Out_Unit,*)'ready is not correct'
   STOP 'sfelevl'
ENDIF

READ(iunit,9,END=2)nstate
IF(ieng==1)THEN
   nchanl = nstate
ENDIF

READ(iunit,*,END=2)(evalues(ienergy),ienergy=istart, istart+nstate-1)
READ(iunit,4,END=2)cont
!IF(cont/='CONTINUE')THEN
!   WRITE(Out_Unit,*)'cont is not correct'
!   STOP 'sfelevl'
!ENDIF

ieng = ieng + 1
IF(neng>mxeng)THEN
   WRITE(Out_Unit,*)'neng is greater than mxeng: neng,mxeng = ', neng, mxeng
   STOP 'sfelevl'
ENDIF

istart = istart + nchanl
IF(istart+nchanl>mxeng*mxstate)THEN
   WRITE(Out_Unit,*)'istart+nchanl>mxstate*neng in SFLevels:',istart,nchanl, mxstate, neng
   STOP 'sfelevl'
ENDIF

GO TO 30
2 neng = ieng - 1

OPEN(Unit=ounit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'GraphicsOut/wgdatAdiabatic.rbw',form='formatted',status='unknown')
   CALL wgdat(engvals, evalues, nchanl, neng, ounit)
CLOSE(ounit)

OPEN(Unit=ounit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'GraphicsOut/pvdatAdiabatic.txt',form='formatted',status='unknown')
   CALL pvdat(engvals, evalues, nchanl, neng, ounit)
CLOSE(ounit)

CALL Date_And_Time(today, hour, curzone, dtvalues)
CALL Print_Date_And_Time(Today, Hour, CurZone, DtValues, Out_Unit)

CLOSE(iunit)
CLOSE(ounit)

DEALLOCATE(engvals, evalues)
WRITE(Out_Unit,*)'Leaving:', ProcName

CLOSE(Out_Unit)
OPEN(Unit=Out_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'Output/BasisForScatt.txt',Form='FORMATTED',ACCESS='APPEND')
WRITE(Out_Unit,*)'Completed SFLevels'
WRITE(Out_Unit,*)

RETURN
8 FORMAT(2x,a5)
9 FORMAT(i6)
3 FORMAT(1x,f10.5)
5 FORMAT(1x,e14.7)
4 FORMAT(2x,a8)
END
