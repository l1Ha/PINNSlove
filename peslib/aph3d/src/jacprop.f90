SUBROUTINE jacprop (njacobi, rmat, eye11, w, vecnew, tstore, vecold, xsq, &
                    enddelve, finish)
USE Numeric_Kinds_Module
!
!------------------------------------------------------------------
!  integrate the aph hyperspherical region. This routine integrates
!  the coupled-channel equations from start to endaph.
!
!  On entering:
!      njacobi           number of coupled-channels in the Delves
!                       hyperspherical region.  Used for the APH
!                       to Delves transformation only.
!      rmat             Wigner R-matrix which is the ratio of the
!                       wavefunction to its derivative.  Initially
!                       this matrix is a (njacobi by njacobi)
!                       real matrix.
!      eye11            work matrix (njacobi by njacobi)
!      w                work matrix (njacobi by njacobi)
!      vecnew           work matrix (njacobi by njacobi)
!      tstore           work matrix (njacobi by njacobi)
!      vecold           work matrix (njacobi by njacobi)
!      xsq              work array  (njacobi)
!      enddelve         starting distance.
!      finish           ending distance.
!
!  On exit:
!      rmat             Wigner R-matrix in the Jacobis basis. This
!                       matrix is a (njacobi by njacobi) real matrix.
!--------------------------------------------------------------------
USE FileUnits_Module
IMPLICIT NONE
REAL(Kind=WP_Kind) enddelve, finish
LOGICAL little, medium, full
INTEGER njacobi, ithcall, ithsub
REAL(Kind=WP_Kind) rmat(*), eye11(njacobi,njacobi)
REAL(Kind=WP_Kind) w(njacobi,njacobi)
REAL(Kind=WP_Kind) tstore(njacobi,njacobi), vecold(njacobi,njacobi)
REAL(Kind=WP_Kind) vecnew(njacobi,njacobi), xsq(njacobi)
EXTERNAL popt, mxout, hybridph
!------------------------------------------------------------------
! Determine printing options.
!------------------------------------------------------------------
DATA little/.false./,medium/.false./,full/.false./, ithcall/0/,ithsub/0/
CALL popt('jacprop ',little,medium,full,ithcall,ithsub)
IF(full)THEN
!------------------------------------------------------------------
! Print the incoming R-matrix if desired.
!------------------------------------------------------------------
    WRITE(Out_Unit,*)' R-matrix at the ', 'Jacobi of the Jacobi region, ', 'enddelve = ', enddelve
    CALL MxOut(rmat, njacobi, njacobi)
ENDIF
!----------------------------------------------------------------------
! integrate in jacobi coordinates.
!----------------------------------------------------------------------
CALL hybridph (njacobi, rmat, eye11, w, vecnew, tstore, vecold, xsq, enddelve, finish)
WRITE(*,*)"Completed hybridph"
WRITE(Out_Unit,*)"Completed hybridph"
!------------------------------------------------------------------
! Print the calculated R-matrix if desired.
!------------------------------------------------------------------
IF(full)THEN
    WRITE(Out_Unit,*)' R-matrix at the ', 'Jacobi of the Jacobi region, ', 'finish = ', finish
    CALL MxOut(rmat, njacobi, njacobi)
ENDIF
RETURN
END
