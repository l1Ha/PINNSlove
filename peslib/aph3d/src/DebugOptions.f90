SUBROUTINE DebugOptions(Out_Unit)
USE Prntng_Module
USE SubList_Module
INTEGER IPR, IT, Out_Unit, PrintTimes
INTEGER, PARAMETER:: PrintMax=1
SAVE PrintTimes
DATA PrintTimes/0/

IF(PrintTimes<PrintMax)THEN
   WRITE(Out_Unit,*)
   WRITE(Out_Unit,*)"Start: DebugOptions"
   WRITE(Out_Unit,*)"OPTION  =", OPTION
   WRITE(Out_Unit,*)"IPRTALL =", IPRTALL
   WRITE(Out_Unit,*)"NSUBS   =", NSUBS
   DO IPR=1,NSUBS
      WRITE(Out_Unit,'(A,I3,A,A,A,I3,A,I3,A,3I3)')"SUBS(",IPR,")=",SUBS(IPR)," IPRT(",IPR,")=",IPRT(IPR)," ITER=",(ITER(IT,IPR),IT=1,3)
   ENDDO
   WRITE(Out_Unit,*)"End: DebugOptions"
ENDIF
PrintTimes=PrintTimes+1

RETURN
END SUBROUTINE DebugOptions