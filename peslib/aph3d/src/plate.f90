SUBROUTINE plate (nvertex, vertice, nrsfs, poly_p, filename)
USE FileUnits_Module
!----------------------------------------------------------------------------------
!
!                   LIGHTWAVE OBJECT: PLATE (plate_)
!
! Runs up to 99 different colored plates*.  The filename must end in a two
! digit number, which corresponds to the color of the object.  The variables
! are the same as in "lwo_" except
!
!   flname1[] == filename string ("_p.lwo" is added automatically. "keith.sub"
!                is changed to "keith_p.lwo" and "keith" is changed to
!                "k_p.lwo")
!
!   * The object does not have to be a plate, but this program was originally
!     written to give a plate with a different color and surface name than
!     the orignal object.  It could be useful if two obects need to be seperate
!     but veiwed together, perhaps as two different layers in Modeler.
!------------------------------------------------------------------------------

IMPLICIT NONE
CHARACTER(LEN=*) filename
INTEGER nvertex, nrsfs
INTEGER(Kind=IW_Kind) poly_p(nrsfs)
REAL(KIND=WP_Kind) vertice(nvertex)
RETURN
END
