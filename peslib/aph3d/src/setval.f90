SUBROUTINE setval(neig,nrho,track)
USE FileUnits_Module
IMPLICIT NONE
INTEGER i, jrho, neig, nrho, track(neig,NRho)

WRITE(Out_Unit,*)
WRITE(Out_Unit,*)'Initializing array'
DO i=1,neig
   track(i,1)=i
   DO jrho=2,nrho
      track(i,jrho)=0
   ENDDO
ENDDO
WRITE(Out_Unit,*)'Completed initializing array'
WRITE(Out_Unit,*)

RETURN
ENDSUBROUTINE setval
