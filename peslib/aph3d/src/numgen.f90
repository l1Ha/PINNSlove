SUBROUTINE numgen
!
!numgen's purpose -- fills up COMMON/vector/ with numbers
!

!this routine is called by: AphDel_Old or AphDel_New
!this routine calls
!
!-----------------------------------------------------------------------
USE FileUnits_Module
USE vector_Module

IMPLICIT NONE

INTEGER i

DO  i=1,64
igr(i)=i-1
ENDDO
fpn=igr
RETURN
END
