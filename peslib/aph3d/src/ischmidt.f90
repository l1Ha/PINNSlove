SUBROUTINE ischmidt (vect, nrow, ncol, orthog)
USE FileUnits_Module
IMPLICIT NONE
INTEGER :: nrow, ncol, jcol, irow


INTEGER :: ddot, dotp, vect (nrow, ncol), orthog (nrow)
! This routine makes the vector orthog orthogonal to all
! of the vectors in vect.
DO jcol = 1, ncol
dotp = 0
ddot = 0
DO irow = 1, nrow
dotp = dotp + vect (irow, jcol) * orthog (irow)
ddot = ddot + vect (irow, jcol) * vect (irow, jcol)
ENDDO
IF(dotp/=0)THEN
   DO irow = 1, nrow
   orthog (irow) = ddot * orthog (irow) - dotp * vect (irow, jcol)
   ENDDO
ENDIF

ENDDO
RETURN
ENDSUBROUTINE ischmidt
