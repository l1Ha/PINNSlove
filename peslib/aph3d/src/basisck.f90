SUBROUTINE basisck (sovrlap, check, naph)
USE Numeric_Kinds_Module
USE FileUnits_Module
USE Numbers_Module
!
!   $RCSfile: basisck.f,v $ $Revision: 1.15 $
!   $Date: 89/11/16 11:02:51 $
!   $State: Stable $
!

!            P U R P O S E   O F    S U B R O U T I N E

! Determines whether the new basis completely expands the old basis.
! The closer check is to unity the better the completeness.

!            I N P U T    A R G U M E N T S

!  sovrlap   Overlap matrix.
!  naph      Number of surface functions

!            O U T P U T    A R G U M E N T S

!  check     Completeness check. This array will be used to sector
!            sizes.

IMPLICIT NONE
INTEGER iaph, kaph, naph, ithcall, ithsub
REAL(Kind=WP_Kind) sovrlap, check, sum
DIMENSION sovrlap(naph, naph), check(naph)
!-----------------------------------------------------------------------
! Determine printing options.
!-----------------------------------------------------------------------
LOGICAL little, medium, full
DATA ithcall /0/, ithsub /0/
DATA little /.false./, medium /.false./, full /.false./
CALL popt ('basisck ', little, medium, full, ithcall, ithsub)
Little=.True.
!-----------------------------------------------------------------------
! Check the basis the reverse direction.
!-----------------------------------------------------------------------
DO iaph = 1, naph
   sum = Zero
   DO kaph = 1, naph
      sum = sum + sovrlap(kaph, iaph)**2
   ENDDO
   check(iaph) = sum
ENDDO
!-----------------------------------------------------------------------
! WRITE out the completeness check.
!-----------------------------------------------------------------------
IF(little)THEN
   WRITE(Out_unit,*)'Column sum completeness check for the new basis'
   WRITE(Out_unit,*)check
ENDIF
!-----------------------------------------------------------------------
! Check the basis in the forward direction.
!-----------------------------------------------------------------------
DO iaph = 1, naph
   sum = Zero
   DO kaph = 1, naph
      sum = sum + sovrlap(iaph, kaph)**2
   ENDDO
   check(iaph) = sum
ENDDO
!-----------------------------------------------------------------------
! WRITE out the completeness check.
!-----------------------------------------------------------------------
IF(little)THEN
   WRITE(Out_unit,*)'Row sum completeness check for the new basis'
   WRITE(Out_unit,*)check
ENDIF

RETURN
ENDSUBROUTINE basisck
