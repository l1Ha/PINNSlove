SUBROUTINE AsymptoticTest
! Author: Gregory A. Parker, Department of Physics and Astronomy, University of Oklahoma
USE Numeric_Kinds_Module
USE FileUnits_Asymptotic_Module
USE InputFile_Module
USE TotalEng_Module
USE Energy_Module
USE EDeriv_Module
IMPLICIT NONE
LOGICAL KNonSym                                                 ! True is K-Matrix is Non-Symmetric
INTEGER :: OutUnit1=948, OutUnit2=949
INTEGER kenergy
INTEGER, PARAMETER:: NOpen=261                                    ! Number of coupled-states.
CHARACTER(LEN=21), PARAMETER:: ProcName='AsymptoticTest'
CHARACTER(LEN=6) Print_Flag
REAL(KIND=WP_Kind)K_Matx(NOpen,NOpen), KE_Matx(NOpen,NOpen)     ! K-Matrix and Energy derivative of K-Matrix

CALL PoptAsy(ProcName, Print_Flag)
Asymptotic_Test=.True.
OPEN(Unit=Out_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'Output/AsymptoticTest.txt',Form='FORMATTED')
WRITE(Out_Unit,*)'Called    AsymptoticTest'
OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//'Input/AllBasis.txt',form='formatted',status='old')
   REWIND(In_Unit)
   READ(In_Unit,NML=IFileName)
   WRITE(Out_Unit,NML=IFilename)
CLOSE(In_Unit)
WRITE(Out_Unit,'("InputFile=",A)')InputFile(1:LEN(InputFile))

OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile,status='old')

nenergy=76
KNonSym=.False.           ! Generate and use non-symmetric K-Matrix
EDeriv=.true.             ! Generate and use the energy derivative of the K-Matrix
OPEN(Unit=OutUnit1,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/K_Matrix.bin',Form='unformatted', Status='unknown')
OPEN(Unit=OutUnit2,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/KE_Matrix.bin',Form='unformatted', Status='unknown')
!DO kenergy=1,nenergy         ! Test with generated K-Matrix and KE-Matrix
!   Etot=.001d0*kenergy+1.d0-.06d0
!   WRITE(Out_Unit,*)'Generated K_Matrix with KNonSym=', KNonSym, ' and EDeriv=', EDeriv
!   CALL K_GEN_TEST(NOpen, K_Matx, KE_Matx, EDeriv, KNonSym, Etot, kenergy, OutUnit1, OutUnit2)        ! Generate dummy K-Matrix  
!ENDDO
CLOSE(OutUnit1)
CLOSE(OutUnit2)   
CALL Asymptotic_Analysis(NOpen)  !

CLOSE(Out_Unit)
OPEN(Unit=Out_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'Output/BasisForScatt.txt',Form='FORMATTED',ACCESS='APPEND')
WRITE(Out_Unit,*)'Completed AsymptoticTest'
WRITE(Out_Unit,*)

STOP
END SUBROUTINE AsymptoticTest
