SUBROUTINE AnalysisOfCurves
USE Numeric_Kinds_Module
USE FileUnits_Module
USE FileUnits_Common_Module
IMPLICIT NONE
LOGICAL there
LOGICAL, PARAMETER:: dbug=.False.
CHARACTER(LEN=1), PARAMETER :: Jobz='E', Uplo='U'
CHARACTER(LEN=80) SubTitle
INTEGER NCurves, NSectors, i, j, lwork, Info
REAL(Kind=WP_Kind),  ALLOCATABLE:: Ovlp(:,:)     ! Ovlp(NCurves,NCurves)
REAL(Kind=WP_Kind),  ALLOCATABLE:: Elevels(:,:)  ! Elevels(NCurves,READ)
REAL(Kind=WP_Kind),  ALLOCATABLE:: Curve(:,:)    ! Curve(NCurves,NSectors)
REAL(Kind=WP_Kind),  ALLOCATABLE:: RhoVals(:)    ! RhoVals(NSectors)

OPEN(Unit=9837,File=OutDIR(1:LEN(TRIM(OutDIR)))//'Output/AnasisOfCurvesDiab.csv',Form='FORMATTED')  
OPEN(Unit=DiaAdia_Bin_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/DiabaticLevels.bin',form='unformatted')
  READ(DiaAdia_Bin_Unit)Subtitle
  READ(DiaAdia_Bin_Unit)NSectors
    ALLOCATE(RhoVals(NSectors))           ! Allocate RhoVals
  READ(DiaAdia_Bin_Unit)(RhoVals(i),i=1,NSectors)
  READ(DiaAdia_Bin_Unit)NCurves
    ALLOCATE(Ovlp(NCurves,NCurves))       ! Allocate Ovlp
    ALLOCATE(Curve(NCurves,NSectors))     ! Allocate Curve
  READ(DiaAdia_Bin_Unit)((Curve(j,i),i=1,NSectors),j=1,NCurves)
CLOSE(DiaAdia_Bin_Unit)

WRITE(*,*)"Calling PossibleResonances using diabatic curves"
CALL PossibleResonances(NCurves, NSectors, Curve, RhoVals)
CLOSE(9387)

OPEN(Unit=9837,File=OutDIR(1:LEN(TRIM(OutDIR)))//'Output/AnasisOfCurvesAdia.csv',Form='FORMATTED')  
OPEN(Unit=DiaAdia_Bin_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/AdiabaticLevels.bin',form='unformatted')
  READ(DiaAdia_Bin_Unit)Subtitle
  READ(DiaAdia_Bin_Unit)NSectors
  READ(DiaAdia_Bin_Unit)(RhoVals(i),i=1,NSectors)
  READ(DiaAdia_Bin_Unit)NCurves
    ALLOCATE(Elevels(NCurves,NSectors))   ! Allocate Elevels
  READ(DiaAdia_Bin_Unit)((Elevels(j,i),i=1,NSectors),j=1,NCurves)
CLOSE(DiaAdia_Bin_Unit)

WRITE(*,*)"Calling PossibleResonances using adiabatic curves"
CALL PossibleResonances(NCurves, NSectors, Elevels, RhoVals)
CLOSE(9387)

DEALLOCATE(Ovlp)
DEALLOCATE(Elevels)
DEALLOCATE(Curve)
DEALLOCATE(RhoVals)

RETURN
ENDSUBROUTINE AnalysisOfCurves