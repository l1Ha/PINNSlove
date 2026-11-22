SUBROUTINE delprop (ndelve, rmat, eye11, w, vecnew, tstore, vecold, xsq, endaph, enddelve)
!
!------------------------------------------------------------------
!  integrate the aph hyperspherical region. This routine integrates
!  the coupled-channel equations from start to endaph.
!
!  On entering:
!      ndelve           number of coupled-channels in the Delves
!                       hyperspherical region.
!      rmat             Wigner R-matrix which is the ratio of the
!                       wavefunction to its derivative.  Initially
!                       this matrix is a (ndelve by ndelve) real matrix.
!      eye11            work matrix (ndelve by ndelve)
!      w                work matrix (ndelve by ndelve)
!      vecnew           work matrix (ndelve by ndelve)
!      tstore           work matrix (ndelve by ndelve)
!      vecold           work matrix (ndelve by ndelve)
!      xsq              work array  (ndelve)
!      endaph           starting hyperradius.
!      enddelve         ending hyperradius.
!
!  On exit:
!      rmat             Wigner R-matrix in the Delves basis. This
!                       matrix is a (ndelves by ndelves) real matrix.
!--------------------------------------------------------------------
USE FileUnits_Module
USE Numeric_Kinds_Module
USE Boundary_Module
USE BasisRegions_Module
IMPLICIT NONE
REAL(Kind=WP_Kind) endaph, enddelve
LOGICAL little, medium, full
INTEGER ndelve, ithcall, ithsub
! PORT (2026-10-02): the Delves propagation rebuilds Basis_Dist/
! Sector_Boundaries (see Delves_Sector_Grid); the per-energy propagator
! (propph1) still expects the APH-region grid, so save and restore it here.
INTEGER :: nsectors_keep
REAL(Kind=WP_Kind), ALLOCATABLE, SAVE :: bd_keep(:), sb_keep(:)
REAL(Kind=WP_Kind) rmat(*), eye11(ndelve,ndelve)
REAL(Kind=WP_Kind) w(ndelve,ndelve), vecnew(ndelve,ndelve)
REAL(Kind=WP_Kind) tstore(ndelve,ndelve), vecold(ndelve,ndelve), xsq(ndelve)
EXTERNAL popt, mxout, hybridph
!------------------------------------------------------------------
! Determine printing options.
!------------------------------------------------------------------
DATA little/.false./,medium/.false./,full/.false./, ithcall/0/,ithsub/0/
CALL popt('delprop ',little,medium,full,ithcall,ithsub)
!------------------------------------------------------------------
! Print the incoming R-matrix if desired.
!------------------------------------------------------------------
IF(full)THEN
    WRITE(Out_Unit,*)' R-matrix at the start of ', 'the Delves region, ', &
            'endaph = ', endaph
    CALL MxOut(rmat, ndelve, ndelve)
ENDIF
!---------------------------------------------------------------------
! integrate in delves hyperspherical coordinates
!---------------------------------------------------------------------
IF(.NOT.ALLOCATED(bd_keep))THEN
   ALLOCATE(bd_keep(Max_Sector), sb_keep(Max_Sector+1))
ENDIF
nsectors_keep = MIN(NSectors, Max_Sector-1)
bd_keep(1:nsectors_keep)        = Basis_Dist(1:nsectors_keep)
sb_keep(1:nsectors_keep+1)      = Sector_Boundaries(1:nsectors_keep+1)
CALL hybridph (ndelve, rmat, eye11, w, vecnew, tstore, vecold, xsq, endaph, enddelve)
! restore the APH-region sector arrays for the per-energy propagator
Basis_Dist(1:nsectors_keep)        = bd_keep(1:nsectors_keep)
Sector_Boundaries(1:nsectors_keep+1) = sb_keep(1:nsectors_keep+1)
NSectors = nsectors_keep
!------------------------------------------------------------------
! Print the calculated R-matrix if desired.
!------------------------------------------------------------------
IF(full)THEN
    WRITE(Out_Unit,*)' R-matrix at the end of the ','Delves region, ', 'enddelve = ', enddelve
    CALL MxOut(rmat, ndelve, ndelve)
ENDIF
RETURN
!----------------------***end-delprop---------------------------------
END
