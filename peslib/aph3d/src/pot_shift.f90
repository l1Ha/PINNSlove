SUBROUTINE pot_shift (n, rmid, w)
!---------------------------------------------------------------------
! written by g. a. parker
! calculates the interaction matrix
!---------------------------------------------------------------------
USE FileUnits_MODULE
IMPLICIT NONE
LOGICAL little, medium, full
INTEGER stdout, n, ithcll, ithsub, ii, iii
!#include <Inc.unitno>
REAL(Kind=WP_Kind) xktot, rmid, w(n*(n+1)/2)
COMMON /utotal/ xktot

DATA ithcll/0/,ithsub/0/
DATA little/.false./, medium/.false./, full/.false./
CALL popt('pot_shift', little, medium, full, ithcll, ithsub)

ii=0
DO iii=1,n
    ii=ii+iii
    w(ii)=w(ii)-xktot
ENDDO
IF(medium)THEN
    WRITE(Out_Unit,*) 'The potential at r=', rmid
    CALL MxOutD(w, n, n, 1)
ENDIF

RETURN
END
