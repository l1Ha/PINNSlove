SUBROUTINE SortEnergy(AsymptEnergies, EminLoc, n)
USE Numeric_Kinds_Module
USE FileUnits_Asymptotic_Module
USE QAsy_Numbers_Module, ONLY : QNumbr
IMPLICIT NONE
LOGICAL MaskArray(n)
INTEGER(KIND=IW_Kind) :: n, i, EminLoc(n), EnergyIndex(1)
CHARACTER(LEN=10), PARAMETER:: ProcName='SortEnergy'
CHARACTER(LEN=6) Print_Flag
REAL(Kind=WP_Kind) :: AsymptEnergies(n)

CALL PoptAsy(ProcName, Print_Flag)
MaskArray=.True.
!-----------------------------------------------------------------------
! Find Energy Ordering
!-----------------------------------------------------------------------
DO i=1,n
   EnergyIndex=MINLOC(AsymptEnergies,Mask=MaskArray)
   EminLoc(i)=EnergyIndex(1)
   MaskArray(EminLoc(i))=.False.
ENDDO

RETURN
ENDSUBROUTINE SortEnergy