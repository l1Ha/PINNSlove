SUBROUTINE K_Gen_Test(NOpen, K_Matx, KE_Matx, EDeriv, KNonSym, Energy, kenergy, outunit1, outunit2)
USE numeric_kinds_module     ! Use numerical precision module
USE FileUnits_Asymptotic_Module         ! Use fileunits module
USE matrixlabels_module      ! Use matrixlabels module for printing quantum labels
USE QAsy_Numbers_Module
USE TotalEng_Module
IMPLICIT NONE
CHARACTER(LEN=21), PARAMETER:: ProcName='K_Gen_Test'
CHARACTER(LEN=6) Print_Flag

!
! Author: Gregory A. Parker, Department of Physics and Astronomy, University of Oklahoma
!
! Generate a dummy K-Matrix to test the Asymptotic codes.
! This routine is used for test purposes only.
!
! Required Input <=====
!     NOpen          Size of ALL matrices (NOpen by NOpen)
!     EDeriv=True    If ernergy derivatives of the K-Matrix exist
!     KNonSym=True   If s non-symmetric K-Matrix is to be generated.
!
! On return  =====>
!     K_Matx         K-Matrix
!     KE_Matx        Energy derivative of the K-Matrix (If EDeriv=.True.)
! 
! This routine is called by:
!      Asymptotic
! This routine calls:
!      Matrix_Out
LOGICAL, INTENT(IN):: EDeriv                                             ! True if Energy derivatives are calculated.
INTEGER, INTENT(IN):: NOpen                                              ! Number of coupled-states.
LOGICAL, INTENT(IN):: KNonSym                                            ! Generate Non-Symmetric K-Matrix
INTEGER IState, JState, index                                            ! Temporary integer variables
INTEGER, PARAMETER:: JTot=1                                              ! Total Angular Momentum (For this test only)
INTEGER:: kenergy, outunit1, outunit2
CHARACTER(LEN=20) TypeOfMat
REAL (KIND=WP_Kind) x, y                                                 ! Temporary variables
REAL (KIND=WP_Kind) Energy                                               ! Total Energy (For this test only)
REAL (KIND=WP_Kind), INTENT(OUT):: K_Matx(NOpen,NOpen), KE_Matx(NOpen,NOpen)  ! K-Matrix and Energy derivative of K-Matrix


CALL PoptAsy(ProcName, Print_Flag)
WRITE(Out_Unit,*)
WRITE(Out_Unit,*)
WRITE(Out_Unit,*)
WRITE(Out_Unit,*)'Entering:', ProcName
ALLOCATE(QNumbr(Nquant,NOpen))                ! chanl, elec, nvib, jrot, lorb
ALLOCATE(xksq(NOpen))              ! two times the reduced mass times the quantity of the total energy minus the int
ALLOCATE(AsymptEnergies(NOpen))          !
ALLOCATE(wvec(NOpen))  

AsymptEnergies=0.d0

DO IState=1,nopen
   xksq(IState)=Etot-AsymptEnergies(IState)
   wvec(IState)=SQRT(ABS(xksq(IState)))
   QNumbr(1,IState)=MOD(IState,3)+1   !chanl(IState)
   QNumbr(2,IState)=MOD(IState,2)     !elec(IState)
   QNumbr(3,IState)=MOD(IState,4)     !nvib(IState)
   QNumbr(4,IState)=MOD(IState,3)     !jrot(IState)
   QNumbr(5,IState)=MOD(IState,5)     !lorb(IState)
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
! Generate Analytic K-Matrix to test asymptotic analysis
!
DO IState=1,NOpen
   DO JState=1,NOpen
      x=IState+JState
      K_Matx(IState,JState)=TAN(Energy*x)
   ENDDO
ENDDO
CALL Matrix_Out(K_Matx,NOpen,NOpen,'K_Matx','Symmetric K-Matrix', Print_Flag)

!
! Add some asymmetry for testing purposes only.
!
IF(KNonSym)THEN
   DO IState=1,NOpen
      DO JState=1,NOpen
         y=IState-JState
         K_Matx(IState,JState)=K_Matx(IState,JState)+0.01d0*SIN(y)
      ENDDO
   ENDDO
   CALL Matrix_Out(K_Matx,NOpen,NOpen,'K_Matx','Non-symmetric K-Matrix', Print_Flag)
ENDIF

!
! Generate Energy derivative of the Analytic K-Matrix to test asymptotic analysis
!
IF(EDeriv)THEN
   DO IState=1,NOpen
      DO JState=1,NOpen
	     x=IState+JState
         KE_Matx(IState,JState)=x/COS(Energy*x)**2
      ENDDO
   ENDDO
   CALL Matrix_Out(KE_Matx,NOpen,NOpen,'KE_Matx','Symmetric d(K-Matrix)/dE', Print_Flag)
   !
   ! Add some asymmetry for testing purposes only.
   !
   IF(KNonSym)THEN
      DO IState=1,NOpen
         DO JState=1,NOpen
            y=IState-JState
            K_Matx(IState,JState)=K_Matx(IState,JState)+0.01d0*SIN(y)
         ENDDO
      ENDDO
      CALL Matrix_Out(KE_Matx,NOpen,NOpen,'KE_Matx','Non-symmetric d(K-Matrix)/dE', Print_Flag)
    ENDIF
ENDIF

TypeOfMat='K_Matrix'
CALL Write_Matrix2(K_Matx, Nopen, Energy, kenergy, TypeOfMat, OutUnit1)
WRITE(Out_Unit,*)'K_Matrix Open-Open Block'
CALL Matrix_Out (K_Matx, nopen, nopen, 'K_Matx', 'K-Matrix', Print_Flag)

IF(EDeriv)THEN
   TypeOfMat='KE_Matrix'
   CALL Write_Matrix2(KE_Matx, Nopen, Energy, kenergy, TypeOfMat, outunit2)
   WRITE(Out_Unit,*)'Energy Derivative of the K_Matrix'
   CALL Matrix_Out (KE_Matx, nopen, nopen, 'KE_Matx', 'Energy Derivative of the K-Matrix', Print_Flag)
ENDIF

DEALLOCATE(QNumbr)   ! elec, chanl, nvib, jrot, lorb
DEALLOCATE(xksq)     ! two times the reduced mass times the quantity of the total energy minus the int
DEALLOCATE(AsymptEnergies) ! Internal energies 
DEALLOCATE(wvec)     ! wavevector

WRITE(Out_Unit,*)'Leaving:', ProcName
RETURN
END SUBROUTINE K_Gen_Test