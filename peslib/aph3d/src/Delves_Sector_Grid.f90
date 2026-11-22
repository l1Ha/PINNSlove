SUBROUTINE Delves_Sector_Grid(rbegin, rmax)
!-----------------------------------------------------------------------
! PORT (2026-10-02): build a sector grid for the Delves propagation region
! [rbegin, rmax] and store it in the shared arrays Basis_Dist /
! Sector_Boundaries / NSectors (Boundary_Module / BasisRegions_Module).
!
! WHY: the shipped APH3D fills those arrays from the ABM/APH input grid
! (RhoStart .. RhoStartDelves) only.  hybridph's sector loop uses them for
! every non-Jacobi region, so a non-zero Delves region (endaph < enddelve)
! ran off the end of Basis_Dist (SIGSEGV at hybridph:129).  Here the grid is
! rebuilt for the Delves interval with logarithmic spacing (as in the
! reference paper) when the Delves propagation starts.
!
!   rbegin    = endaph    (start of the Delves region, bohr)
!   rmax      = enddelve  (end of the Delves region, bohr)
!   nsect_in  = requested number of sectors (<=0 -> automatic ~0.25 bohr)
!-----------------------------------------------------------------------
USE Numeric_Kinds_Module
USE Boundary_Module
USE BasisRegions_Module
IMPLICIT NONE
REAL(Kind=WP_Kind), INTENT(IN) :: rbegin, rmax
INTEGER :: i, n
REAL(Kind=WP_Kind) :: rlo, rhi, ratio

rlo = MIN(rbegin, rmax)
rhi = MAX(rbegin, rmax)
IF (rhi - rlo < 1.d-10) THEN
   n = 2
ELSE
   n = INT((rhi - rlo)/0.25d0)          ! ~0.25 bohr average spacing
ENDIF
n = MAX(4, MIN(n, Max_Sector - 2))

NSectors = n
IF (rhi - rlo < 1.d-10) THEN
   DO i = 1, n
      Basis_Dist(i) = rlo
   ENDDO
   DO i = 1, n + 1
      Sector_Boundaries(i) = rlo
   ENDDO
ELSE
   ratio = (rhi/rlo)**(1.d0/REAL(n, WP_Kind))
   ! sector edges (log spaced) and centres
   Sector_Boundaries(1) = rlo
   DO i = 2, n + 1
      Sector_Boundaries(i) = Sector_Boundaries(i-1)*ratio
   ENDDO
   Sector_Boundaries(n+1) = rhi
   DO i = 1, n
      Basis_Dist(i) = SQRT(Sector_Boundaries(i)*Sector_Boundaries(i+1))
   ENDDO
ENDIF

RETURN
END SUBROUTINE Delves_Sector_Grid
