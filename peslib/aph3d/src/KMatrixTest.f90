SUBROUTINE KMatrixTest
! Author: Gregory A. Parker, Department of Physics and Astronomy, University of Oklahoma
USE Numeric_Kinds_Module
USE FileUnits_Asymptotic_Module
USE InputFile_Module
USE TotalEng_Module
USE Energy_Module
USE EDeriv_Module
USE QAsy_Numbers_Module
IMPLICIT NONE
LOGICAL KNonSym                                                 ! True is K-Matrix is Non-Symmetric
INTEGER :: OutUnit1=948, OutUnit2=949
INTEGER kenergy, nopen
INTEGER, PARAMETER :: NOpenx=90                                    ! Number of coupled-states.
CHARACTER(LEN=21), PARAMETER:: ProcName='KMatrixTest'
CHARACTER(LEN=6) Print_Flag
REAL(KIND=WP_Kind)URU_Matx(NOpenx,NOpenx), UREU_Matx(NOpenx,NOpenx)     ! K-Matrix and Energy derivative of K-Matrix

NOpen=NOpenx
OPEN(Unit=Out_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'Output/KMatrixTest.txt',Form='FORMATTED')
WRITE(Out_Unit,*)'Called    KMatrixTest'
OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//'Input/AllBasis.txt',form='formatted',status='old')
   REWIND(In_Unit)
   READ(In_Unit,NML=IFileName)
   WRITE(Out_Unit,NML=IFilename)
CLOSE(In_Unit)
WRITE(Out_Unit,'("InputFile=",A)')InputFile(1:LEN(InputFile))

OPEN(Unit=In_Unit,File=TRIM(InputDIR)//TRIM(InputFile),status='old')

nenergy=121
KNonSym=.True.               ! Generate and use non-symmetric K-Matrix
EDeriv=.true.             ! Generate and use the energy derivative of the K-Matrix
OPEN(Unit=OutUnit1,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/URU_Matrix.bin',Form='unformatted', Status='unknown')
OPEN(Unit=OutUnit2,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/UREU_Matrix.bin',Form='unformatted', Status='unknown')
DO kenergy=1,nenergy         ! Test with generated K-Matrix and KE-Matrix
   Etot=.001d0*kenergy+1.d0-.06d0
   WRITE(Out_Unit,*)'Generated K_Matrix with KNonSym=', KNonSym, ' and EDeriv=', EDeriv
   CALL URU_GEN_TEST(NOpen, URU_Matx, UREU_Matx, EDeriv, KNonSym, Etot, kenergy, OutUnit1, OutUnit2)        ! Generate dummy K-Matrix  
ENDDO
DEALLOCATE(Conserved)
DEALLOCATE(ConLabels)
DEALLOCATE(QLabels)
CLOSE(OutUnit1)
CLOSE(OutUnit2)   
CALL KMatrx_Main(NOpen)  !

CLOSE(Out_Unit)
OPEN(Unit=Out_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'Output/BasisForScatt.txt',Form='FORMATTED',ACCESS='APPEND')
WRITE(Out_Unit,*)'Completed KMatrixTest'
WRITE(Out_Unit,*)

STOP
END SUBROUTINE KMatrixTest
