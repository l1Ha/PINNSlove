SUBROUTINE Asymptotic(NOpen, K_Matx, KE_Matx, EDeriv, KNonSym)
USE Numeric_Kinds_Module
USE FileUnits_Asymptotic_Module
USE matrixlabels_Module      ! Use matrixlabels module for printing quantum labels
USE ReadK_Module
USE QAsy_Numbers_Module, ONLY : NQuant, QNumbr
IMPLICIT NONE
CHARACTER(LEN=10), PARAMETER:: ProcName='Asymptotic'
CHARACTER(LEN=6) Print_Flag
!
! Author: Gregory A. Parker, Department of Physics and Astronomy, University of Oklahoma
!
! Calculates S-Matrix, T-Matrix, Probabilities and Eigenphases from the K-Matrix
! IF(EDeriv=.True.) Time-Delays and Lifetimes will be calculated also.
!                    Note: This requires energy derivative of the K-Matrix (KE_Matx).
!
! Required Input <=====
!     NOpen       Size of ALL matrices (NOpen by NOpen)
!     EDeriv=T    If ernergy derivatives of the K-Matrix exist
!     KNonSym=T   If non-symmetric K-Matrix was generated and unitarity is an issue
!     K_Matx      K-Matrix
!     KE_Matx     Energy derivative of the K-Matrix
!
! This routine is called by:
!      Main
! This routine calls:
!      K_Gen_Test
!      S_Matrix
!      Probabilities
!      Unitarity
!      Symmetrize
!      Eigenphase
!      Lifetimes           If EDeriv=True
!      Time_Delays         If EDeriv=True
!
CHARACTER(LEN=5)  CurZone
CHARACTER(LEN=8)  Today
CHARACTER(LEN=10) Hour
CHARACTER(LEN=20) TypeOfMat
LOGICAL, INTENT(IN):: EDeriv                         ! True if Energy derivatives are calculated.
LOGICAL, INTENT(IN):: KNonSym                        ! True if K_Matrix is not symmetric and Non-Unitary S-Matrix is desired
INTEGER, INTENT(IN):: NOpen                          ! Number of coupled-states.
INTEGER I,J
INTEGER  DtValues(8)
REAL(Kind=WP_Kind), INTENT(INOUT):: K_Matx(NOpen,NOpen) ! K-Matrix
REAL(Kind=WP_Kind), INTENT(INOUT):: KE_Matx(NOpen,NOpen)! Energy derivative of K-Matrix
!
! Allocatable storage
!
REAL(Kind=WP_Kind), ALLOCATABLE:: S_Real(:,:), S_Imag(:,:)          ! S-Matrix (real and imaginary parts)
REAL(Kind=WP_Kind), ALLOCATABLE:: T_Real(:,:), T_Imag(:,:)          ! T-Matrix (real and imaginary parts)
REAL(Kind=WP_Kind), ALLOCATABLE:: Prob(:,:)                         ! Probabilities
REAL(Kind=WP_Kind), ALLOCATABLE:: Phase(:), Phase_Tran(:,:)         ! Eigenphases and Transformation Matrix

! These arrarys are only needed if EDeriv=.true.
REAL(Kind=WP_Kind), ALLOCATABLE:: SE_Real(:,:), SE_Imag(:,:)        ! Energy derivative of the S-Matrix
REAL(Kind=WP_Kind), ALLOCATABLE:: TE_Real(:,:), TE_Imag(:,:)        ! Energy derivative of the T-Matrix
REAL(Kind=WP_Kind), ALLOCATABLE:: ProbE(:,:)                        ! Energy derivative of Probabilities
REAL(Kind=WP_Kind), ALLOCATABLE:: PhaseE(:), PhaseE_Tran(:,:)       ! Energy derivative of the Eigenphases and Transformation Matrix
REAL(Kind=WP_Kind), ALLOCATABLE:: Q_Real(:,:), Q_Imag(:,:)          ! Lifetimes (real and imaginary parts)
REAL(Kind=WP_Kind), ALLOCATABLE:: TD_Real(:,:), TD_Imag(:,:)        ! Time-Delay Matrix
!
! Open output files
!
OPEN(Unit=eigunit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'GraphicsOut/Eigenphases.rbw')
IF(EDeriv)THEN
   OPEN(Unit=derunit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'GraphicsOut/EigenphasesE.rbw')
ENDIF
CALL Date_And_Time(today, hour, curzone, dtvalues)
CALL Print_Date_And_Time(Today, Hour, CurZone, DtValues, Out_Unit)
WRITE(Out_Unit,*)
WRITE(Out_Unit,*)
WRITE(Out_Unit,*)
WRITE(Out_Unit,*)'Entering:', ProcName
WRITE(Out_Unit,*)'Energy Derivatives=',EDeriv
WRITE(Out_Unit,*)'Number of coupled-states=',NOpen
LabelRows=.True.
LabelColumns=.True.
CALL PoptAsy(ProcName, Print_Flag)


!
! Allocate storage for Asymptotic matrices.
!
ALLOCATE(S_Real(NOpen,NOpen),S_Imag(NOpen,NOpen))      ! S-Matrix
ALLOCATE(T_Real(NOpen,NOpen),T_Imag(NOpen,NOpen))      ! T-Matrix
ALLOCATE(Prob(NOpen,NOpen))                            ! Probability Matrix
ALLOCATE(Phase(NOpen),Phase_Tran(NOpen,NOpen))         ! Eigenphases and Eigenphase Transformation Matrix

!
! Allocate storage for Asymptotic matrices that require energy derivatives.
!
IF(EDeriv)THEN
   ALLOCATE(SE_Real(NOpen,NOpen),SE_Imag(NOpen,NOpen)) ! S-Matrix
   ALLOCATE(TE_Real(NOpen,NOpen),TE_Imag(NOpen,NOpen)) ! T-Matrix
   ALLOCATE(ProbE(NOpen,NOpen))                        ! Probability Matrix
   ALLOCATE(PhaseE(NOpen),PhaseE_Tran(NOpen,NOpen))    ! Eigenphases and Eigenphase Transformation Matrix
   ALLOCATE(Q_Real(NOpen,NOpen),Q_Imag(NOpen,NOpen))   ! Lifetime Matrix
   ALLOCATE(TD_Real(NOpen,NOpen),TD_Imag(NOpen,NOpen)) ! Time-Delay Matrix
ENDIF


!
! Determine S-Matrix and Probabilities from unsymmeterized K-Matrix.
!
IF(KNonSym)THEN
   CALL S_Matrix(NOpen,K_Matx,S_Real,S_Imag,KE_Matx,SE_Real,SE_Imag,EDeriv)             ! Get S-Matrix
   CALL Probabilities(NOpen,S_Real,S_Imag,Prob,SE_Real,SE_Imag,ProbE,EDeriv)            ! Get Probabilities
   CALL Unitarity(NOpen, S_Real, S_Imag, SE_Real, SE_Imag, EDeriv)                      ! Get Unitarity
ENDIF

!
! Symmeterize the K-Matrix
! Determine S-Matrix, T-Matrix, Probabilities and Eigenphases from symmeterized K-Matrix.
!
CALL Symmetrize(NOpen, K_Matx, KE_Matx, EDeriv)
CALL S_Matrix(NOpen,K_Matx,S_Real,S_Imag,KE_Matx,SE_Real,SE_Imag,EDeriv)                ! Get S-Matrix
CALL T_Matrix(NOpen,S_Real,S_Imag,SE_Real,SE_Imag,T_Real,T_Imag,TE_Real,TE_Imag,EDeriv) ! Get T-Matrix
CALL Probabilities(NOpen,S_Real,S_Imag,Prob,SE_Real,SE_Imag,ProbE,EDeriv)               ! Get Probabilities
CALL Unitarity(NOpen, S_Real, S_Imag, SE_Real, SE_Imag, EDeriv)                         ! Get Unitarity
CALL Eigenphase(NOpen,K_Matx,Phase,Phase_Tran,KE_Matx,PhaseE,PhaseE_Tran,EDeriv)        ! Get Eigenphases

! PORT (2026-10-02): dump the complex S-matrix as text so that elastic
! amplitudes |1-S|^2 (needed for the Wigner threshold analysis) can be
! post-processed; proba.csv only contains |S|^2.
OPEN(Unit=87, File=OutDIR(1:LEN(TRIM(OutDIR)))//'Output/Smatrix.txt', &
     Status='UNKNOWN', Position='APPEND')
WRITE(87,'("# energy(Ha)=",ES24.16," Nopen=",I6)') etot, nopen
! PORT 2026-10-10: the channel ordering of this dump is the K-matrix
! open-channel order (ascending original channel index).  Write the
! quantum numbers (chanl, elect, nvib, jrot, lorb) of each column so the
! post-processing can identify channels without guessing the ordering.
WRITE(87,'("# QN chanl elect nvib jrot lorb")')
DO j=1,nopen
   WRITE(87,'("# QN",5I6)') (QNumbr(i,j), i=1,5)
ENDDO
DO j=1,nopen
   DO i=1,nopen
      WRITE(87,'(2ES24.16)') S_Real(i,j), S_Imag(i,j)
   ENDDO
ENDDO
CLOSE(Unit=87)

!
! Determine Eigenphases, Lifetimes and Time-Delays from symmeterized K-Matrix.
!
IF(EDeriv)THEN
   CALL Lifetimes(NOpen,S_Real,S_Imag,SE_Real,SE_Imag,Q_Imag,Q_Real)                   ! Get Lifetimes
   CALL Time_Delays(NOpen,S_Real,S_Imag,SE_Real,SE_Imag,TD_Imag,TD_Real)               ! Get Time-Delays
ENDIF


TypeOfMat='KSym_Matrix'
CALL Write_Matrix2(K_Matx, nopen, etot, ithmat, TypeOfMat, KSym_Unit)
TypeOfMat='SReal_Matrix'
CALL Write_Matrix2(S_Real, nopen, etot, ithmat, TypeOfMat, SReal_Unit)
TypeOfMat='SImag_Matrix'
CALL Write_Matrix2(S_Imag, nopen, etot, ithmat, TypeOfMat, SImag_Unit)
TypeOfMat='Prob_Matrix'
!WRITE(832,'(2x,i5,e14.8,2x,i5,A)')Nopen, etot, ithmat, typeofmat
!Write(832,'(11000(E12.4,","))')etot*autoev,((Prob(i,j),i=1,nopen),j=1,nopen)
CALL Write_Matrix2(Prob, nopen, etot, ithmat, TypeOfMat, Probil_Unit)
TypeOfMat='Phase_Matrix'
CALL Write_Vector2(Phase, Nopen, etot, ithmat, TypeOfMat, Phase_Unit)
TypeOfMat='PhaseTran_Matrix'
CALL Write_Matrix2(Phase_Tran, nopen, etot, ithmat, TypeOfMat, PhaseTran_Unit)
                   
                   
IF(Ederiv)THEN
   TypeOfMat='KESym_Matrix'
   CALL Write_Matrix2(KE_Matx, nopen, etot, ithmat, TypeOfMat, KESym_Unit)
   TypeOfMat='SEReal_Matrix'
   CALL Write_Matrix2(SE_Real, nopen, etot, ithmat, TypeOfMat, SEReal_Unit)
   TypeOfMat='SEImag_Matrix'
   CALL Write_Matrix2(SE_Imag, nopen, etot, ithmat, TypeOfMat, SEImag_Unit)
   TypeOfMat='ProbE_Matrix'
   CALL Write_Matrix2(ProbE, nopen, etot, ithmat, TypeOfMat, ProbilE_Unit)
   TypeOfMat='PhaseE_Matrix'
   CALL Write_Vector2(PhaseE, Nopen, etot, ithmat, TypeOfMat, PhaseE_Unit)
   TypeOfMat='PhaseTranE_Matrix'
   CALL Write_Matrix2(PhaseE_Tran, nopen, etot, ithmat, TypeOfMat, PhaseTranE_Unit)
   TypeOfMat='QReal_Matrix'
   CALL Write_Matrix2(Q_Real, nopen, etot, ithmat, TypeOfMat, QReal_Unit)
   TypeOfMat='QImag_Matrix'
   CALL Write_Matrix2(Q_Imag, nopen, etot, ithmat, TypeOfMat, QImag_Unit)
   TypeOfMat='TDelayReal_Matrix'
   CALL Write_Matrix2(TD_Real, nopen, etot, ithmat, TypeOfMat, TDelayReal_Unit)
   TypeOfMat='TDelayImag_Matrix'
   CALL Write_Matrix2(TD_Imag, nopen, etot, ithmat, TypeOfMat, TDelayImag_Unit)
ENDIF                 


WRITE(Out_Unit,*)'Deallocate Temporary Storage in asymptotic'

!
! Deallocate storage for Asymptotic matrices.
!
DEALLOCATE(Prob)               ! Probability Matrix
DEALLOCATE(S_Real,S_Imag)      ! S-Matrix
DEALLOCATE(T_Real,T_Imag)      ! T-Matrix
DEALLOCATE(Phase,Phase_Tran)   ! Eigenphases and Eigenphase Transformation Matrix

!
! Deallocate storage for Asymptotic matrices that require energy derivatives.
!
IF(EDeriv)THEN
   DEALLOCATE(ProbE)                 ! Probability Matrix
   DEALLOCATE(SE_Real,SE_Imag)       ! S-Matrix
   DEALLOCATE(TE_Real,TE_Imag)       ! T-Matrix
   DEALLOCATE(PhaseE,PhaseE_Tran)    ! Eigenphases and Eigenphase Transformation Matrix
   DEALLOCATE(Q_Real,Q_Imag)         ! Lifetime Matrix
   DEALLOCATE(TD_Real,TD_Imag)       ! Time-Delay Matrix
ENDIF

CALL Date_And_Time(today, hour, curzone, dtvalues)
CALL Print_Date_And_Time(Today, Hour, CurZone, DtValues, Out_Unit)
WRITE(Out_Unit,*)
WRITE(Out_Unit,*)
WRITE(Out_Unit,*)'Leaving:', ProcName

RETURN
ENDSUBROUTINE Asymptotic
