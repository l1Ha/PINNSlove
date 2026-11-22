SUBROUTINE check_int (nxval, ndaf, nchanl, mval)
USE FileUnits_Module
!-----------------------------------------------------------------------
! This routine was written by G. A. Parker
! If you find an error or have an improvement please send a messge to
! Parker@ou.edu
!-----------------------------------------------------------------------
IMPLICIT NONE
LOGICAL :: data_ok


INTEGER :: nxval, ndaf, nchanl, mval
!-----------------------------------------------------------------------
! check all input integers to make sure that they are reasonable.
!-----------------------------------------------------------------------

data_ok = .true.
IF(nxval<1)THEN
   data_ok = .false.
ENDIF
IF(ndaf<0)THEN
   data_ok = .false.
ENDIF
IF(nchanl<1)THEN
   data_ok = .false.
ENDIF
IF(mval<2.or. (mval / 2) * 2/=mval)THEN
   WRITE(Out_Unit, * ) 'mval must be an even INTEGER > 2'
   data_ok = .false.
ENDIF
IF(data_ok)THEN
   WRITE(Out_Unit, * ) 'Input integers are reasonable'
   WRITE(Out_Unit, * ) 'nxval= ', nxval
   WRITE(Out_Unit,  * ) 'ndaf=  ', ndaf
   WRITE(Out_Unit, * ) 'nchanl=', nchanl
   WRITE(Out_Unit,  * ) 'mval=  ', mval
ELSE
   WRITE(Out_Unit, * ) 'Input integers are NOT valid'
   WRITE(Out_Unit, * ) 'nxval, ndaf, nchanl, mval'
   WRITE(Out_Unit, * ) nxval, ndaf, nchanl, mval
   STOP 'check_int'

ENDIF
RETURN
ENDSUBROUTINE check_int
