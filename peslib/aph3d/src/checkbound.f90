SUBROUTINE checkbound (xvals, nxval, delta)
USE Numeric_Kinds_Module
USE FileUnits_Module
!-----------------------------------------------------------------------
! This routine was written by G. A. Parker
! If you find an error or have an improvement please send a messge to
! Parker@ou.edu
!-----------------------------------------------------------------------
IMPLICIT NONE
LOGICAL :: data_ok
INTEGER :: nxval, ix
REAL(Kind=WP_Kind) :: delta, xvals (nxval)
REAL(Kind=WP_Kind), PARAMETER :: tol = 1.d-6
!-----------------------------------------------------------------------
! check to make sure that the input data is uniform.
!-----------------------------------------------------------------------
data_ok = .true.
delta = xvals (2) - xvals (1)
WRITE(Out_Unit, * ) 'Deltax=', delta
WRITE(Msg_Unit, * ) 'Deltax=', delta
DO ix = 2, nxval
   IF(ABS(delta-(xvals(ix)-xvals(ix-1))) >tol*ABS(delta))THEN
      WRITE(Out_Unit, * ) 'Error: Grid not uniform'
      WRITE(Out_Unit, * ) 'ix=', ix, delta , (xvals (ix) - xvals (ix - 1) ) 
      WRITE(Msg_Unit, * ) 'Error: Grid not uniform'
      WRITE(Msg_Unit, * ) 'ix=', ix, delta , (xvals (ix) - xvals (ix - 1) ) 
      WRITE(Msg_Unit, * ) 'ix=', ix,ABS(delta-(xvals(ix)-xvals(ix-1))),tol*ABS(delta)
      data_ok = .false.
   ENDIF
 ENDDO
IF(.NOT.data_ok)THEN
   WRITE(Out_Unit, * ) 'Error: Grid not uniform'
   STOP 'CheckBound'
ELSE
   WRITE(Out_Unit, * ) 'input xvals OK DATA is uniform'
ENDIF

RETURN
ENDSUBROUTINE checkbound
