SUBROUTINE aphprop (naph, rmat, eye11, w,  &
             vecnew, tstore, vecold, xsq, start, endaph, &
             ndelve,xk,nidm)
!
!------------------------------------------------------------------
!  integrate the aph hyperspherical region. This routine integrates
!  the coupled-channel equations from start to endaph.
!
!  On entering:
!      naph             number of coupled-channels in the APH
!                       hyperspherical region.
!      ndelve           number of coupled-channels in the Delves
!                       hyperspherical region.  Used for the APH
!                       to Delves transformation only.
!      rmat             Wigner R-matrix which is the ratio of the
!                       wavefunction to its derivative.  Initially
!                       this matrix is a (naph by naph) real matrix.
!      eye11            work matrix (naph by naph)
!      w                work matrix (naph by naph)
!      vecnew           work matrix (naph by naph)
!      tstore           work matrix (naph by naph)
!      vecold           work matrix (naph by naph)
!      xsq              work array  (naph)
!      start            starting hyperradius (bohr)
!      endaph           ending hyperradius (bohr)
!      jrmax            value of the largest rotational quantum
!                       number used in the Delves region.
!  On exit:
!      rmat             Wigner R-matrix in the Delves basis. This
!                       matrix is a (ndelves by ndelves) real matrix.
!--------------------------------------------------------------------
USE FileUnits_Module
USE rotmax_Module
IMPLICIT NONE
!#include <Inc.unitno>
REAL(Kind=WP_Kind) rmat, eye11, w, vecnew, tstore
REAL(Kind=WP_Kind) vecold, xsq, start, endaph, xk
LOGICAL little, medium, full
INTEGER naph, ndelve, ithcall, ithsub, nidm
!#include <Inc.rotmax>
DIMENSION rmat(*), eye11(naph,naph), w(naph,naph), vecnew(naph,naph), &
    tstore(naph,naph), vecold(naph,naph), xsq(naph), xk(*)
LOGICAL inext
EXTERNAL popt, mxout, aphget1, hybridph, upsiln
!------------------------------------------------------------------
! Determine printing options.
!------------------------------------------------------------------
DATA little/.false./,medium/.false./,full/.false./, ithcall/0/,ithsub/0/
CALL popt('aphprop ',little,medium,full,ithcall,ithsub)
!--------------------------------------------------------------------
! Get the quantum numbers associated with each channel.
!--------------------------------------------------------------------
inext=.false.
CALL aphget1(naph)
!--------------------------------------------------------------------
IF(full)THEN
    WRITE(Out_Unit,*)' R-matrix at the start of the APH region, ', 'start = ', start
    CALL MxOut(rmat, naph, naph)
ENDIF
!--------------------------------------------------------------------
! Propagate the equations from start to endaph.
! At the completion of this CALL
!    rmat   contains the Wigner R-matrix in the APH basis
!           evaluated at the distance endaph.
!--------------------------------------------------------------------
IF(start<=endaph)THEN
   CALL hybridph (naph, rmat, eye11, w, vecnew, tstore, vecold, xsq, start, endaph)
ENDIF
!----------------------------------------------------------------------
!  Print the r-matrix in the APH basis if desired.
!----------------------------------------------------------------------
IF(full)THEN
   WRITE(Out_Unit,*)'  R-matrix after APH propagation, endaph = ',  endaph
   CALL MxOut(rmat, naph, naph)
ENDIF

RETURN
!---------------------***end-aphprop-----------------------------------
ENDSUBROUTINE aphprop
