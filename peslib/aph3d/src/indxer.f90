SUBROUTINE indxer(index, ms, ithrow, ithcol, nrow, ncol)
USE FileUnits_Module
!
!  p u r p o s e   o f    s u b r o u t i n e

!---------------------------------------------------------------------
!  this routine finds the indexes of a matrix for printing purposes.
!  this routine called from routine mxoutl.
!
!  on entering
!  ms       storage option parameter.
!  ms=0 matrix stored as a general matrix.
!  ms=1 matrix stored as a symmetric matrix.
!  ithrow   row index.
!  ithcol   initial column index.
!  nrow     number of rows in the matrix
!  ncol     number of columns of the matrix
!
!  on exit
!  index    contains the indexes for 10 positions of the matrix.
!----------------------------------------------------------------------

!  i n p u t    a r g u m e n t s
!
!  index
!  ms
!  ithrow
!  ithcol
!  nrow
!  ncol

USE FileUnits_Module
IMPLICIT NONE
INTEGER index, ithrow, ithcol, ind, jthcol, ms, nrow, ncol
DIMENSION index(1:10)
!----------------------------------------------------------------------
!  the matrix is stored as a general matrix.
!----------------------------------------------------------------------
IF(ms == 0)THEN
   DO 10 ind=1, 10
      index(ind)=ithrow+(ithcol+ind-2)*nrow
 10      CONTINUE
!----------------------------------------------------------------------
!  the matrix is stored as a symmetric matrix.
!----------------------------------------------------------------------
ELSEIF(ms == 1)THEN
   DO 20 ind=1, 10
      jthcol=ithcol+ind-1
      IF(jthcol <= ithrow)THEN
         index(ind)=ind+ithrow*(ithrow-1)/2
      ELSE
         index(ind)=index(ind-1)+jthcol-1
      ENDIF
 20      CONTINUE
ELSE
   WRITE(Out_Unit,*)'execution stopped in routine indexer.'
   WRITE(Out_Unit,*)'the storage parameter is neither ', '0 nor 1: ms=', ms
   STOP 'indxer'
ENDIF
END
