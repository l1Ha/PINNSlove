SUBROUTINE Tracker(neig, nrho, ovlp, track, rhovals)
USE FileUnits_Module
USE Numeric_Kinds_Module
IMPLICIT NONE
INTEGER neig, nrho, jrho, track(neig,NRho)
REAL(Kind=WP_Kind) rho, rholast
REAL(Kind=WP_Kind) ovlp(neig,neig), rhovals(nrho)

WRITE(Out_Unit,*)
WRITE(Out_Unit,*)'Reading overlaps'
DO jrho=2,nrho
   IF(jrho==nrho)WRITE(Out_Unit,*)'jrho=',jrho
   CALL ReadOvlp(neig,ovlp,rholast,rho)
   IF(jrho==2)THEN
      rhovals(1)=rholast
      WRITE(Out_Unit,*)'jrho=',jrho-1,rholast
   ENDIF
   rhovals(jrho)=rho
   CALL connect(neig,nrho,jrho,ovlp,track)
ENDDO

WRITE(Out_Unit,*)'jrho=',nrho,rho
WRITE(Out_Unit,*)'Completed reading overlaps'
WRITE(Out_Unit,*)

RETURN
ENDSUBROUTINE tracker
