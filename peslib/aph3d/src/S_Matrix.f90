SUBROUTINE S_Matrix(NOpen, K_Matx, S_Real, S_Imag, KE_Matx, SE_Real, SE_Imag, EDeriv)
USE Numeric_Kinds_Module
USE FileUnits_Asymptotic_Module
IMPLICIT NONE
CHARACTER(LEN=8), PARAMETER:: ProcName='S_Matrix'
CHARACTER(LEN=6) Print_Flag

! Author: Gregory A. Parker, Department of Physics and Astronomy, University of Oklahoma
!
! Calculate S-Matrix from the K-Matrix
!
! Required Input <=====
!     NOpen       Size of ALL matrices (NOpen by NOpen)
!     EDeriv=True If ernergy derivatives of the K-Matrix exist
!     K_Matx      K-Matrix
!     KE_Matx     Energy derivative of the K-Matrix
!
! On return  =====>
!     S_Real      Real part of the S-Matrix
!     S_Imag      Imaginary part of the S-Matrix
!     SE_Real     Energy derivative of the real part of the S-Matrix
!     SE_Imag     Energy derivative of the imaginary part of the S-Matrix
!
! This routine is called by:
!      Asymptotic
! This routine calls:
!      Matrix_Out
!      DGESV     a LAPACK routine for solving a system of linear equations
INTEGER, INTENT(IN):: NOpen
INTEGER  Info
INTEGER, ALLOCATABLE:: ipiv(:)
LOGICAL, INTENT(IN):: EDeriv
REAL(Kind=WP_Kind), INTENT(IN)::  K_Matx(NOpen,NOpen),  KE_Matx(NOpen,NOpen)
REAL(Kind=WP_Kind), INTENT(OUT):: S_Real(NOpen,NOpen),  S_Imag(NOpen,NOpen)
REAL(Kind=WP_Kind), INTENT(OUT):: SE_Real(NOpen,NOpen), SE_Imag(NOpen,NOpen)
REAL(Kind=WP_Kind), ALLOCATABLE:: Temp(:,:), UnitMat(:,:)

WRITE(Out_Unit,*)
WRITE(Out_Unit,*)
WRITE(Out_Unit,*)
WRITE(Out_Unit,*)'Entering:', ProcName
CALL PoptAsy(ProcName, Print_Flag)
ALLOCATE (Temp(NOpen,NOpen), UnitMat(NOpen,NOpen), ipiv(NOpen))

S_Real=0.d0
S_Imag=0.d0

CALL Unit_Matrix(NOpen, UnitMat)
Temp=TRANSPOSE(UnitMat+Matmul(K_Matx,K_Matx))
S_Imag=2.d0*TRANSPOSE(K_Matx)
!-------------------------------------------------------------------
! Temp now contains TRANSPOSE[1+(K_Matx)*(K_Matx)]
! S_Imag now contains 2*TRANSPOSE(K_Matx)
!-------------------------------------------------------------------
CALL dgesv(NOpen, NOpen, Temp, NOpen, ipiv, S_Imag, NOpen, Info)
S_Imag=TRANSPOSE(S_Imag)
IF(Info/=0)THEN
   WRITE(Msg_Unit,*)'Cannot Solve S_Matrix Equations'
   WRITE(Out_Unit,*)'Cannot Solve S_matrix Equations'
   STOP 'S_Matrix'
ENDIF
!-------------------------------------------------------------------
! S_Imag now contains the imaginary part of the S_matrix.
!-------------------------------------------------------------------

S_Real=UnitMat-MatMul(S_Imag,K_Matx)
!-------------------------------------------------------------------
! S_Real now contains the real part of the s-matrix
!-------------------------------------------------------------------

IF(EDeriv)THEN
   Temp=UnitMat+MatMul(K_Matx,K_Matx)
   SE_Imag=-MatMul(MatMul(KE_Matx,K_Matx)+MatMul(K_Matx,KE_Matx),S_Imag)+2*KE_Matx
   CALL dgesv(NOpen,NOpen,Temp,NOpen,ipiv,SE_Imag,NOpen,Info)
   IF(Info/=0)THEN
      WRITE(Msg_Unit,*)'Cannot Solve S_Matrix Derivative Equations'
      WRITE(Out_Unit,*)'Cannot Solve S_Matrix Derivative Equations'
      STOP 'S_Matrix EDeriv=.True.'
   ENDIF
   !-------------------------------------------------------------------
   ! SE_Imag now contains the derivative of the imaginary part of the
   ! S_matrix with respect to energy.
   !-------------------------------------------------------------------
   SE_Real=-MatMul(SE_Imag,K_Matx)-MatMul(S_Imag,KE_Matx)
   !-------------------------------------------------------------------
   ! SE_Real now contains the derivative of the real part of the
   ! S_matrix with respect to energy.
   !-------------------------------------------------------------------
ENDIF
!
! Print S-Matrix and if possible the energy derivative of the S-Matrix
!
CALL Matrix_Out(S_Real,NOpen,NOpen,'S_Real','Real Part of S-Matrix', Print_Flag)
CALL Matrix_Out(S_Imag,NOpen,NOpen,'S_Imag','Imaginary Part of S-Matrix', Print_Flag)
IF(EDeriv)THEN
   CALL Matrix_Out(SE_Real,NOpen,NOpen,'SE_Real','d(Real Part of S-Matrix)/dE', Print_Flag)
   CALL Matrix_Out(SE_Imag,NOpen,NOpen,'SE_Imag','d(Imaginary Part of S-Matrix)/dE', Print_Flag)
ENDIF

WRITE(Out_Unit,*)'Deallocate Temporary Storage in S_Matrix'
DEALLOCATE (Temp, UnitMat, ipiv)

WRITE(Out_Unit,*)'Leaving:', ProcName
RETURN
ENDSUBROUTINE S_Matrix
