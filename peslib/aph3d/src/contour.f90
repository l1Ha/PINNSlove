SUBROUTINE contour (nvertex, vert, npolygon, poly, filename, dimensions, ncontours, seperation, cont_start)
!------------------------------------------------------------------------------------
!
!                   LIGHTWAVE OBJECT with CONTOURS (contours_)
!
! This program creates either a 2D or 3D contour map of an object.  The
! variables that differ from "lwo_" are as follows:
!
!   flname1[] == filename string ("_c.lwo" is added automatically. "keith.sub"
!                is changed to "keith_c.lwo" and "keith" is changed to
!                "k_c.lwo")
!
!   D == number of dimensions.  2 correspondes to a flat contour map, while
!        all other numbers correspond to a 3D contour map.
!
!   ncont == the number of contour lines (up to 999, in theory)
!
!   rise == the seperation between each contour line
!
!   start == the initial contour value
!---------------------------------------------------------------------------------------------
USE Numeric_Kinds_Module
IMPLICIT NONE
CHARACTER(LEN=*) filename
INTEGER nvertex, npolygon, ncontours, dimensions, poly(npolygon)
REAL(KIND=WP_Kind) vert(3,nvertex), seperation
REAL(KIND=WP_Kind) cont_start

RETURN
END
