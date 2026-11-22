SUBROUTINE Adiabatic(NSectors)
USE Numeric_Kinds_Module
USE FileUnits_Module
USE APHChl_Module, ONLY       : neigmin_sum, naph
USE RhoSur_Module, Only       : Curve
USE QAsy_Numbers_Module, ONLY : Nquant, Qnumbr, Qlabels, AsymptEnergies
USE Converge_Module, ONLY     : NGood
IMPLICIT NONE
LOGICAL there
LOGICAL, PARAMETER:: dbug=.False.
INTEGER NCurves, NSectors, i, j
INTEGER, ALLOCATABLE:: track(:,:)                !track(NCurves,NSectors)
REAL(Kind=WP_Kind),  ALLOCATABLE:: ovlp(:,:)     !ovlp(NCurves,NCurves)
REAL(Kind=WP_Kind),  ALLOCATABLE:: elevels(:,:)  !elevels(NCurves,NSectors)
!REAL(Kind=WP_Kind),  ALLOCATABLE:: curve(:,:)    !curve(NCurves,NSectors)
REAL(Kind=WP_Kind),  ALLOCATABLE:: rhovals(:)    !rhovals(NSectors)

ncurves=naph !tmpmodgregparker maybe neigmin_sum
ncurves=Min(NGood,naph)  !238
ncurves=ngood  !tmpmodgregparker
IF(.not.ALLOCATED(QNumbr))ALLOCATE(QNumbr(NQuant, ncurves))                ! chanl, elec, nvib, jrot, lorb
IF(.not.ALLOCATED(AsymptEnergies))ALLOCATE(AsymptEnergies(ncurves))        !
IF(.not.ALLOCATED(QLabels))  ALLOCATE(QLabels(NQuant))
!IF(.not.ALLOCATED(ConLabels))ALLOCATE(ConLabels(NConserved))
!IF(.not.ALLOCATED(Conserved))ALLOCATE(Conserved(NConserved))
WRITE(Out_Unit,*)'Called    Adiabatic'
CLOSE(Out_Unit)
OPEN(Unit=Out_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'Output/Adiabatic.txt',Form='FORMATTED')

ALLOCATE(track(NCurves,NSectors))
ALLOCATE(ovlp(NCurves,NCurves))
ALLOCATE(elevels(NCurves,NSectors))
ALLOCATE(curve(NCurves,NSectors))
ALLOCATE(rhovals(NSectors))

INQUIRE(file=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/Ovrlp.bin', exist=there)
IF(THERE)THEN
   OPEN(Unit=DiaAdia_Bin_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/Ovrlp.bin',status='old',form='unformatted')
ELSE
   WRITE(Out_Unit,*)'Error in opening BinOut/Ovrlp_All'
   WRITE(Out_Unit,*)'Stopping in Adiabatic'
   STOP 'Adiabatic'
ENDIF

!------------------------------------------------------------------------
! READ in energy levels and rhovalues.
!------------------------------------------------------------------------
CALL ReadLevels(NCurves, NSectors, elevels, rhovals)

!------------------------------------------------------------------------
! Initialize the tracking array.
!------------------------------------------------------------------------
CALL setval(NCurves, NSectors, track)
IF(Dbug)THEN
   WRITE(Out_Unit,*)'track after setval='
   DO i=1,NCurves
      WRITE(Out_Unit,'(22i3)')(track(i,j),j=1,NSectors)
   ENDDO
ENDIF
!------------------------------------------------------------------------
! Follow the curves.
!------------------------------------------------------------------------
CALL tracker(NCurves,NSectors,ovlp,track,rhovals)
IF(Dbug)THEN
   WRITE(Out_Unit,*)'track after tracker='
   DO i=1,NCurves
      WRITE(Out_Unit,'(22i3)')(track(i,j),j=1,NSectors)
   ENDDO
ENDIF

!------------------------------------------------------------------------
! Determine how many segments are connected
!------------------------------------------------------------------------
CALL CheckAdiabatic(NCurves,NSectors,track)
IF(Dbug)THEN
   WRITE(Out_Unit,*)'track after checkadiabatic='
   DO i=1,NCurves
      WRITE(Out_Unit,'(22i3)')(track(i,j),j=1,NSectors)
   ENDDO
ENDIF

!------------------------------------------------------------------------
! Print out DATA for energy correlation plots
!------------------------------------------------------------------------
CALL PrCorr(NCurves,NSectors,track,elevels,rhovals,curve)
IF(Dbug)THEN
   WRITE(Out_Unit,*)'track PrCorr='
   DO i=1,NCurves
      WRITE(Out_Unit,'(22i3)')(track(i,j),j=1,NSectors)
   ENDDO
ENDIF


CLOSE(Out_Unit)
OPEN(Unit=Out_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'Output/BasisForScatt.txt',Form='FORMATTED',ACCESS='APPEND')
WRITE(Out_Unit,*)'Completed Adiabatic'
WRITE(Out_Unit,*)

RETURN
ENDSUBROUTINE Adiabatic
