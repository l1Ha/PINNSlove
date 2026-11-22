

!-------------------------------------------------------------------
! SUBROUTINE: WRITEtime
!
! package   : CID
!
! Language  : Fortran 90
!
! author    : F. Colavecchia (flavioc@lanl.gov)
!
! date      :              version:
! revision  :              version:
!
! purpose   :  Get the boundary matrix in Jacobi coordinates
!              open basis.
!
! input     :    qnh -> Index of bound states in the basis
!                qsf -> Index of bound states in the SF
!              nchanl-> Number of channels
!
! output    :    a1  -> Regular bessel function matrix at rhomatch
!                b1  -> Irregular bessel function matrix at rhomatch
!
! modules   :
!
!
! common    :
!
!
! notes     :
!
!-------------------------------------------------------------------

SUBROUTINE WRITEtime(procname,time)

! Modules

USE Time_Conv_Module
USE fileunits_Module

IMPLICIT NONE

! Input / Output variables

CHARACTER(LEN=*) procname
CHARACTER(LEN=13), PARAMETER :: tblank='             '
CHARACTER(LEN=1), PARAMETER :: sblank=' '
REAL(Kind=WP_Kind) time
! Internals

WRITE(time_unit,'(a13,a10,a1,f10.3)') tblank,procname,sblank,time

RETURN
END
