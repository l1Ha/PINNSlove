SUBROUTINE graf(start, delta, psi, npts, filename, nds, ider)
USE Numeric_Kinds_Module
USE FileUnits_Module
USE FileUnits_OneDim_Module
!-----------------------------------------------------------------------
! This routine was written by G. A. Parker
! IF you find an error or have an improvement please send a messge to
! Parker@ou.edu 
!-----------------------------------------------------------------------
IMPLICIT NONE
CHARACTER*(*) filename
INTEGER npts, ip, nds, nstart, ider
REAL(KIND=WP_Kind) start, delta, psi(nds:npts)
REAL(KIND=SP_Kind) xval, pltf

!----------------------------------------------------------------------- 
! outputs results for subsequent plotting.
!-----------------------------------------------------------------------

OPEN(Unit=Graf_Unit,file=TRIM(OutDIR)//"PlotVee/"//filename)
IF(nds.eq.1)THEN
   nstart=1
ELSE IF(nds.eq.0)THEN
   nstart=-npts
ELSE
   STOP 'Graf'
ENDIF

DO ip=nstart,npts
   xval=start+(ip-1)*delta
   IF(ip.lt.0)THEN
      pltf=psi(-ip)*(-1)**ider
   ELSE
      pltf=psi(ip)
   ENDIF
   WRITE(Graf_Unit,'(1x,6(1pe14.6,","))')xval,pltf
ENDDO
CLOSE(Unit=Graf_Unit)

RETURN
END SUBROUTINE graf
