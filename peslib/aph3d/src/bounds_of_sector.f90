SUBROUTINE bounds_of_sector (rho_last, rho_now, rho_next, sector_left,  &
                             sector_right, start, finish, isector,      &
                             first_sector, last_sector)
USE boundary_Module
USE BasisRegions_Module
USE fileunits_Module
USE nrho_Module
USE Time_Conv_Module

! I/O variables
IMPLICIT NONE
REAL(Kind=WP_Kind) rho_now, rho_last, rho_next, sector_left, sector_right
REAL(Kind=WP_Kind) start, finish, Rho_Basis
INTEGER isector, first_sector, last_sector, isect

CHARACTER(LEN=21), PARAMETER:: procname='bounds_of_sector'
Logical DEBUG
REAL(Kind=WP_Kind) timeinfo

IF(isector>max_sector)THEN
   WRITE(msg_unit,*) 'isector    =',isector
   WRITE(msg_unit,*) 'max_sector =',max_sector
   CALL errormsg(msg_unit,procname, 'isector gt max_sector')
ENDIF

!-----------------------------------------------------------------------
! On the first CALL determine first_sector and last_sector
!----------------------------------------------------------------------
IF(isector<=0)THEN
   first_sector=1
   last_sector=NSectors
   DO isect=2,NSectors
      IF((basis_dist(isect-1)+basis_dist(isect))/2.d0<=start)THEN
            first_sector=isect
      ENDIF
      IF((basis_dist(isect)+basis_dist(isect+1))/2.d0<finish)THEN
            last_sector=isect
      ENDIF
   ENDDO
   WRITE(Out_Unit,*)'first_sector=',first_sector,' last_sector=',last_sector
   RETURN
ENDIF

rho_basis=basis_dist(isector)
WRITE(Out_Unit,*)'SECTOR=',isector,'  rho_basis=',rho_basis
WRITE(msg_unit,*)'SECTOR=',isector,'  rho_basis=',rho_basis

!-----------------------------------------------------------------------
! Compute the left and right boundaries of the sector.
!-----------------------------------------------------------------------
IF(isector>first_sector.AND.isector<last_sector)THEN
   sector_left=sector_boundaries(isector)
   sector_right=sector_boundaries(isector+1)
   rho_last=basis_dist(isector-1)
   rho_now =basis_dist(isector)
   rho_next=basis_dist(isector+1)
ELSEIF (isector==first_sector.AND.isector<last_sector)THEN
   sector_left=start
   sector_right=sector_boundaries(isector+1)
   rho_last=basis_dist(isector)
   rho_now =basis_dist(isector)
   rho_next=basis_dist(isector+1)
ELSEIF (isector>first_sector.AND.isector==last_sector)THEN
   sector_left=sector_boundaries(isector)
   sector_right=finish
   rho_last=basis_dist(isector-1)
   rho_now =basis_dist(isector)
   rho_next=basis_dist(isector)
ELSEIF (isector==first_sector.AND.isector==last_sector)THEN
   sector_left=start
   sector_right=finish
   rho_last=basis_dist(isector)
   rho_now =basis_dist(isector)
   rho_next=basis_dist(isector)
ENDIF
RETURN
END
