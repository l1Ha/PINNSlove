PROGRAM Program_TestLifetimes
USE Numeric_Kinds_Module       ! Set numeric kinds
USE FileUnits_Asymptotic_Module           ! Fileunit numbers
USE MatrixLabels_Module        ! Use matrixlabels module for printing quantum labels
USE physical_Constants_Module  ! Contains physical constants
   !
   ! This program calculates lifetimes and eigen lifetimes.
   !

IMPLICIT NONE
INTEGER, PARAMETER:: inunit1=81      ! Unit for S_Real            [Real part of the S-matrix           ]
INTEGER, PARAMETER:: inunit2=82      ! Unit for S_Imag            [Imaginary part of the S-matrix      ]
INTEGER, PARAMETER:: inunit3=83      ! Unit for S_Real_deriv      [Real part of the (dS-matrix/dE)     ]
INTEGER, PARAMETER:: inunit4=84      ! Unit for S_Imag_deriv      [Imaginary part of the (dS-matrix/dE)]
INTEGER, PARAMETER:: NTest=10
INTEGER NState, I

REAL(KIND=WP_Kind), ALLOCATABLE:: S_Real(:,:), Q_real(:,:)
REAL(KIND=WP_Kind), ALLOCATABLE:: S_Imag(:,:), Q_imag(:,:)
REAL(KIND=WP_Kind), ALLOCATABLE:: SE_Real(:,:), SE_Imag(:,:)

OutDir=TRIM(OutputRoot)//"Testing/Test_Lifetimes/"       ! Output directory 
InputDir=TRIM(InputRoot)//"Testing/Test_Lifetimes/"  ! Input directory

OPEN (UNIT=Out_Unit, file=TRIM(OutDir)//'Output/Lifetimes.txt', STATUS='unknown')
OPEN (UNIT=EigLifePltUnit, file=TRIM(OutDir)//'GraphicsOut/EigLife.csv', STATUS='unknown')
CALL NumericInformation
NState=2
DO I=1,1
   ALLOCATE(S_Real(NState,NState))         ! Real part of the S-Matrix
   ALLOCATE(S_Imag(NState,NState))         ! Imaginary part of the S-Matrix
   ALLOCATE(SE_Real(NState,NState))        ! Energy derivative of the Real(S-Matrix)
   ALLOCATE(SE_Imag(NState,NState))        ! Energy derivative of the Imag(S-Matrix)
   ALLOCATE(Q_real(NState,NState))         ! Real part of the Lifetime matrix (Q)
   ALLOCATE(Q_imag(NState,NState))         ! Imaginary part of the Lifetime matrix (Q)
   CALL LifetimeTestData(NState, S_Real, S_Imag, SE_Real, SE_Imag)
   CALL Lifetimes (NState, S_Real, S_Imag, SE_Real, SE_Imag, Q_real, Q_imag)
   DEALLOCATE(S_Real)
   DEALLOCATE(Q_real)
   DEALLOCATE(S_Imag)
   DEALLOCATE(Q_imag)
   DEALLOCATE(SE_Real)
   DEALLOCATE(SE_Imag)
END DO

WRITE(Msg_Unit,*)"Successful Calculation of Lifetimes"
WRITE(Out_Unit,*)"Successful Calculation of Lifetimes"
STOP "Program_LifetimeTest Completed Successfully"
END PROGRAM Program_TestLifetimes