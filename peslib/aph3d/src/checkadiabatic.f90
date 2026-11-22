SUBROUTINE CheckAdiabatic(neig,nrho,track)
USE FileUnits_Module
IMPLICIT NONE
INTEGER maxpos,npos,i,jrho,neig,nrho,track(neig,nrho)

maxpos=neig*nrho
npos=0
DO i=1,neig
   DO jrho=1,nrho
      IF(track(i,jrho)/=0)npos=npos+1
   ENDDO
ENDDO

WRITE(Out_Unit,*)'npos=',npos,'   maxpos=',maxpos

RETURN
ENDSUBROUTINE CheckAdiabatic
