SUBROUTINE Combine_SFLevels(SFunType, NType)
USE FileUnits_Module
USE Boundary_Module
IMPLICIT NONE
LOGICAL there
INTEGER ieng, neng, mxeng, ienergy, nstate, iunit, istart, nchanl, IType
INTEGER NType
INTEGER, PARAMETER:: mxstate=1000   ! PORT 2026-10-04: raised 600->1000
CHARACTER(LEN=21), PARAMETER:: ProcName='Combine_SFLevels'
CHARACTER(LEN=5 ) readx, ready
CHARACTER(LEN=8 ) cont
CHARACTER(LEN=30) In_File
CHARACTER(LEN=*) SFunType(NType)
REAL, ALLOCATABLE:: engvals(:), evalues(:)
EXTERNAL wgdat
mxeng=Max_Sector+5
ALLOCATE(EngVals(mxeng), EValues(mxeng*mxstate))

ieng = 1
istart = 1
iunit=21

OPEN(Unit=SFLevel_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'GraphicsOut/Sfelevl_All.rbw',Form='formatted',Status='Unknown')
DO 2 IType=1,NType
   In_File=TRIM('GraphicsOut/Sfelevl_'//SFunType(IType))//'.rbw'
   INQUIRE (File=OutDIR(1:LEN(TRIM(OutDIR)))//In_File, exist=there)
   IF(there)THEN
      OPEN(Unit=iunit,File=OutDIR(1:LEN(TRIM(OutDIR)))//In_File,status='old')
   ELSE
      WRITE(Out_Unit,*)In_File, 'No sfelevl file exist: Stopping in Combine_SFLevels'
      WRITE(Out_Unit,*)OutDIR(1:LEN(TRIM(OutDIR)))//In_File
      WRITE(Msg_Unit,*)In_File, 'No sfelevl file exist: Stopping in Combine_SFLevels'
      WRITE(Msg_Unit,*)OutDIR(1:LEN(TRIM(OutDIR)))//In_File
      STOP 'Combine_SFLevels'
   ENDIF

   DO While(IEng<=MxEng)
      READ(IUnit,8,END=2)
      READ(iunit,8,END=2)readx
      IF(readx/='readx')THEN
         WRITE(Out_Unit,*)'readx is not correct:',readx
         WRITE(Out_Unit,*)'Stopping in All_Combine'
         STOP 'Combine_SFLevels'
      ENDIF
      WRITE(SFLevel_Unit,8)readx

      READ(iunit,9,END=2)neng
      IF(neng/=1)THEN
         WRITE(Out_Unit,*)'neng/=1; neng=',neng
         STOP 'Combine_SFLevels'
      ENDIF
      WRITE(SFLevel_Unit,9)neng

      READ(iunit,3,END=2)engvals(ieng)
      READ(iunit,8,END=2)ready
      IF(ready/='ready')THEN
         WRITE(Out_Unit,*)'ready is not correct'
         STOP 'Combine_SFLevels'
      ENDIF
      WRITE(SFLevel_Unit,3)engvals(ieng)
      WRITE(SFLevel_Unit,8)ready

      READ(iunit,9,END=2)nstate
      IF(ieng==1)THEN
         nchanl = nstate
      ENDIF
      WRITE(SFLevel_Unit,9)nstate

      READ(iunit,*,END=2)(evalues(ienergy),ienergy=istart, istart+nstate-1)
      READ(iunit,4,END=2)cont
      WRITE(SFLevel_Unit,*)(evalues(ienergy),ienergy=istart, istart+nstate-1)
      WRITE(SFLevel_Unit,4)cont

      ieng = ieng + 1
      IF(neng>mxeng)THEN
         WRITE(Out_Unit,*)'neng is greater than mxeng: neng,mxeng = ', neng, mxeng
         STOP 'Combine_SFLevels'
      ENDIF

      istart = istart + nchanl
      IF(istart+nchanl>mxeng*mxstate)THEN
        WRITE(Out_Unit,*)'istart+nchanl>mxstate*neng in Combine_Sfevels:',istart,nchanl, mxstate, mxeng
        STOP 'Combine_SFLevels'
      ENDIF
   ENDDO
   CLOSE(iunit)
2 ENDDO
CLOSE(SFLevel_Unit)
DEALLOCATE(EngVals, EValues)

RETURN

8 FORMAT(2x,a5)
9 FORMAT(i6)
3 FORMAT(1x,f10.5)
5 FORMAT(1x,e14.7)
4 FORMAT(2x,a8)

ENDSUBROUTINE Combine_SFLevels
