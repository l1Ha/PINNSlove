SUBROUTINE Graf_DAF(start, delta, psi, npts, filename, nds, ider)
USE FileUnits_Module
USE Numeric_Kinds_Module
!-----------------------------------------------------------------------
! This routine was written by G. A. Parker
! If you find an error or have an improvement please send a messge to
! Parker@ou.edu
!-----------------------------------------------------------------------
IMPLICIT NONE
CHARACTER(LEN=*) filename
INTEGER npts, ip, nds, nstart, ider, ithcall
REAL(Kind=WP_Kind) start, delta, psi(nds:npts)
REAL(KIND=WP_Kind) xval, pltf
!data ithcall/0/
!ithcall=ithcall+1

!-----------------------------------------------------------------------
! outputs results for subsequent plotting.
!-----------------------------------------------------------------------

OPEN(Unit=Graf_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//filename)
!WRITE(*,*)"Ithcall=",ithcall, "filename=", TRIM(filename)
IF(nds==1)THEN
   nstart=1
ELSEIF(nds==0)THEN
   nstart=-npts
ELSE
   STOP 'Graf'
ENDIF

DO ip=nstart,npts
   xval=start+(ip-1)*delta
   IF(ip<0)THEN
      pltf=psi(-ip)*(-1)**ider
   ELSE
      pltf=psi(ip)
   ENDIF
   WRITE(Graf_Unit,'(1x,6(1pe14.6,","))')xval,pltf
ENDDO
CLOSE(Unit=Graf_Unit)

RETURN
ENDSUBROUTINE Graf_DAF
