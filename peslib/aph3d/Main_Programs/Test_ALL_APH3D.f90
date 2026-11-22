PROGRAM Test_ALL_APH3D
INTEGER :: CSTAT, ESTAT
CHARACTER(100) :: CMSG
INTEGER(4) ErrNum, IErrNo
CALL EXECUTE_COMMAND_LINE("DIR", ExitStat=Estat, cmdstat=cstat, cmdmsg=cmsg)

!CALL EXECUTE_COMMAND_LINE("E:\VisualFortran\APH3D\APH3D_PKH3\APH3D_PKH3\x64\Debug\APH3D_PKH3.exe", ExitStat=Estat, cmdstat=cstat, cmdmsg=cmsg)
I=System("E:\VisualFortran\APH3D\APH3D_PKH3\APH3D_PKH3\x64\Debug\APH3D_PKH3.exe")
I=SystemQQ("E:\VisualFortran\APH3D\APH3D_PKH3\APH3D_PKH3\x64\Debug\APH3D_PKH3.exe")
IF (CSTAT > 0) THEN
  WRITE(*,*)"Execution Failed", TRIM(CMSG)
ELSE IF (CSTAT < 0) THEN
  WRITE(*,*)"Command not supported"
ELSE
  WRITE(*,*)"Completed with status=",ESTAT
END IF
   
   
ENDPROGRAM Test_ALL_APH3D