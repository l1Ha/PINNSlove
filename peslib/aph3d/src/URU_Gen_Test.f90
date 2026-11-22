SUBROUTINE URU_Gen_Test(NOpen, URU_Matx, UREU_Matx, EDeriv, KNonSym, Energy, kenergy, outunit1, outunit2)
USE numeric_kinds_module     ! Use numerical precision module
USE FileUnits_Asymptotic_Module         ! Use fileunits module
USE MatrixLabels_module      ! Use matrixlabels module for printing quantum labels
USE QAsy_Numbers_Module
USE TotalEng_Module
!USE Momentum_Module
IMPLICIT NONE
INTEGER, ALLOCATABLE :: IDummy(:)
CHARACTER(LEN=21), PARAMETER:: ProcName='URU_Gen_Test'
CHARACTER(LEN=6) Print_Flag
!
! Author: Gregory A. Parker, Department of Physics and Astronomy, University of Oklahoma
!
! Generate a dummy R-Matrix to test the Asymptotic codes.
! This routine is used for test purposes only.
!
! Required Input <=====
!     NOpen          Size of ALL matrices (NOpen by NOpen)
!     EDeriv=True    If ernergy derivatives of the R-Matrix exist
!     KNonSym=True   If s non-symmetric R-Matrix is to be generated.
!
! On return  =====>
!     URU_Matx         R-Matrix
!     UREU_Matx        Energy derivative of the R-Matrix (If EDeriv=.True.)
! 
! This routine is called by:
!      Asymptotic
! This routine calls:
!      Matrix_Out
LOGICAL, INTENT(IN):: EDeriv                                             ! True if Energy derivatives are calculated.
INTEGER, INTENT(IN):: NOpen                                              ! Number of coupled-states.
LOGICAL, INTENT(IN):: KNonSym                                            ! Generate Non-Symmetric R-Matrix
INTEGER IState, JState, index                                            ! Temporary integer variables
INTEGER, PARAMETER:: JTot=0                                              ! Total Angular Momentum (For this test only)
INTEGER:: kenergy, outunit1, outunit2
INTEGER:: Chanl(100), elec(100), nvib(100), jrot(100), lorb(100)
CHARACTER(LEN=20) TypeOfMat, approx, proj, symm
REAL (KIND=WP_Kind) x, y                                                 ! Temporary variables
REAL (KIND=WP_Kind) Energy                                               ! Total Energy (For this test only)
REAL (KIND=WP_Kind), INTENT(OUT):: URU_Matx(NOpen,NOpen), UREU_Matx(NOpen,NOpen)  ! R-Matrix and Energy derivative of R-Matrix


WRITE(Out_Unit,*)
WRITE(Out_Unit,*)
WRITE(Out_Unit,*)
WRITE(Out_Unit,*)'Entering:', ProcName
ALLOCATE(QNumbr(NQuant,NOpen))                ! chanl, elec, nvib, jrot, lorb
ALLOCATE(xksq(NOpen))              ! two times the reduced mass times the quantity of the total energy minus the int
ALLOCATE(AsymptEnergies(NOpen))          !
ALLOCATE(wvec(NOpen))
ALLOCATE(IDummy(NOpen))

AsymptEnergies=0.d0

DO IState=1,nopen
   xksq(IState)=Etot-AsymptEnergies(IState)
   wvec(IState)=SQRT(ABS(xksq(IState)))
   QNumbr(1,IState)=MOD(IState,2)     !chanl(IState)
   Chanl(IState)=QNumbr(1,IState)
   QNumbr(2,IState)=MOD(IState,3)+1   !elec(IState)
   elec(IState)=QNumbr(1,IState)
   QNumbr(3,IState)=MOD(IState,4)     !nvib(IState)
   nvib(IState)=QNumbr(1,IState)
   QNumbr(4,IState)=MOD(IState,3)     !jrot(IState)
   jrot(IState)=QNumbr(1,IState)
   QNumbr(5,IState)=MOD(IState,5)     !lorb(IState)
   lorb(IState)=QNumbr(1,IState)
ENDDO

!
! Generate Dummy Quantum Labels (for testing only)
!
NQLabels=NQuant
IF(.NOT.ALLOCATED(QStates))THEN
   ALLOCATE(QStates(NQLabels,NOpen))
ENDIF

IF(NQLabels>0)THEN
    RowNames=' '
    ColumnNames='  '
    DO Index=1,NQLabels
       RowNames=TRIM(RowNames)//'  '//TRIM(QLabels(Index)(1:1))
       ColumnNames=TRIM(ColumnNames)//'  '//TRIM(QLabels(Index)(1:1))
    ENDDO
ENDIF

LabelRows=.True.
LabelColumns=.True.
QStates=QNumbr


!
! Generate Analytic R-Matrix to test asymptotic analysis
!
DO IState=1,NOpen
   DO JState=1,NOpen
      x=IState+JState
      URU_Matx(IState,JState)=TAN(Energy*x)
   ENDDO
ENDDO
CALL Matrix_Out(URU_Matx,NOpen,NOpen,'URU_Matx','Symmetric R-Matrix', Print_Flag)

!
! Add some asymmetry for testing purposes only.
!
IF(KNonSym)THEN
   DO IState=1,NOpen
      DO JState=1,NOpen
         y=IState-JState
         URU_Matx(IState,JState)=URU_Matx(IState,JState)+0.01d0*SIN(y)
      ENDDO
   ENDDO
   CALL Matrix_Out(URU_Matx,NOpen,NOpen,'URU_Matx','Non-symmetric R-Matrix', Print_Flag)
ENDIF

!
! Generate Energy derivative of the Analytic R-Matrix to test asymptotic analysis
!
IF(EDeriv)THEN
   DO IState=1,NOpen
      DO JState=1,NOpen
	     x=IState+JState
         UREU_Matx(IState,JState)=x/COS(Energy*x)**2
      ENDDO
   ENDDO
   CALL Matrix_Out(UREU_Matx,NOpen,NOpen,'UREU_Matx','Symmetric d(R-Matrix)/dE', Print_Flag)
   !
   ! Add some asymmetry for testing purposes only.
   !
   IF(KNonSym)THEN
      DO IState=1,NOpen
         DO JState=1,NOpen
            y=IState-JState
            URU_Matx(IState,JState)=URU_Matx(IState,JState)+0.01d0*SIN(y)
         ENDDO
      ENDDO
      CALL Matrix_Out(UREU_Matx,NOpen,NOpen,'UREU_Matx','Non-symmetric d(R-Matrix)/dE', Print_Flag)
    ENDIF
ENDIF
ETOT=Energy
TypeOfMat='URU_Matrix'
!CALL Write_Matrix2(URU_Matx, Nopen, Energy, kenergy, TypeOfMat, OutUnit1)
CALL Write_Matrix2old(URU_Matx, AsymptEnergies, Chanl, elec, nvib, jrot, lorb, nopen, Etot, kenergy, &
                         TypeOfMat, jtot, proj, approx, symm, OutUnit1)
WRITE(Out_Unit,*)'URU_Matrix Open-Open Block'
CALL Matrix_Out (URU_Matx, nopen, nopen, 'URU_Matx', 'R-Matrix', Print_Flag)

IF(EDeriv)THEN
   TypeOfMat='UREU_Matrix'
   !CALL Write_Matrix2(UREU_Matx, Nopen, Energy, kenergy, TypeOfMat, outunit2)
   CALL Write_Matrix2old(UREU_Matx, AsymptEnergies, Chanl, elec, nvib, jrot, lorb, nopen, Etot, kenergy, &
                         TypeOfMat, jtot, proj, approx, symm, OutUnit2)
   WRITE(Out_Unit,*)'Energy Derivative of the URU_Matrix'
   CALL Matrix_Out (UREU_Matx, nopen, nopen, 'UREU_Matx', 'Energy Derivative of the R-Matrix', Print_Flag)
ENDIF

DEALLOCATE(QNumbr)   ! elec, chanl, nvib, jrot, lorb
DEALLOCATE(xksq)     ! two times the reduced mass times the quantity of the total energy minus the int
DEALLOCATE(AsymptEnergies) ! Internal energies 
DEALLOCATE(wvec)     ! wavevector
DEALLOCATE(IDummy)

WRITE(Out_Unit,*)'Leaving:', ProcName
RETURN
END SUBROUTINE URU_Gen_Test