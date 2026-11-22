SUBROUTINE ReadOvlp(neig,ovlp,rholast,rho)
USE Numeric_Kinds_Module
USE FileUnits_Module
IMPLICIT NONE
INTEGER ithrho, neigmin, neig, jbasis, ibasis
REAL(Kind=WP_Kind) ovlp(neig,neig)
REAL(Kind=WP_Kind) rholast, rho

READ(DiaAdia_Bin_Unit) ithrho, rholast, rho
READ(DiaAdia_Bin_Unit) neigmin, rholast, rho

IF(neigmin/=neig)THEN
  WRITE(Out_Unit,*) 'neigmin/=neig',neigmin,neig
  STOP
ENDIF

DO jbasis=1,neigmin
   READ(DiaAdia_Bin_Unit)(ovlp(ibasis,jbasis),ibasis=1,neigmin)  ! OVR_FEM_unit
ENDDO

RETURN
ENDSUBROUTINE ReadOvlp
