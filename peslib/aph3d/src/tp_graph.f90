SUBROUTINE tp_graph(nchanl, nxval, xvals, upot, eig, ntp,filename)
USE FileUnits_Module
USE Numeric_Kinds_Module
USE DiatomicPot_Module, ONLY: DiatomicPot
IMPLICIT NONE
LOGICAL, PARAMETER:: Dbug=.False.
INTEGER nchanl, nxval, ntp(nchanl), i, k
CHARACTER(LEN=*) :: filename
REAL(Kind=WP_Kind) :: xvals (nxval), upot (nxval), eig (nchanl)

REAL(Kind=WP_Kind), ALLOCATABLE :: xtp (:, :)

ALLOCATE (xtp (204, nchanl) )
OPEN(Unit=TP_Graph_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//filename,status='unknown')
WRITE(TP_Graph_Unit,'(3(A13,","))')" ","TPoint ","Eigenvalue"
DO i=1,nchanl
   CALL findtp ( nxval, upot, xvals, eig(i), ntp(i), xtp(1,i))
   IF(Dbug)WRITE(Out_Unit,*)'For i=',i,' eig(i)=',eig(i),'  there are ntp=', ntp(i),' turning points'
   DO k=1,ntp(i)
      IF(Dbug)WRITE(Out_Unit,*)'xtp,eval=',xtp(k,i),eig(i)
   ENDDO
   IF(ntp(i)>0)THEN
      DO k=1,ntp(i),2
         WRITE(TP_Graph_Unit, '(1x,A13,",",2(1pe13.5,","))')" ", xtp (k, i), eig (i)
         IF(k+1<=ntp(i))THEN
            WRITE(TP_Graph_Unit,'(1x,A13,",",2(1pe13.5,","))')" ", xtp (k + 1, i), eig (i)
            IF(k/=ntp(i))WRITE(TP_Graph_Unit,*)'&'
         ENDIF
      ENDDO
   ENDIF
ENDDO
CLOSE(Unit=TP_Graph_Unit)

DEALLOCATE(xtp)

RETURN
ENDSUBROUTINE tp_graph
