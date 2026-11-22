SUBROUTINE Eigenphase (NOpen, K_Matx, Phase, Phase_Tran, KE_Matx, PhaseE, PhaseE_Tran, EDeriv)
USE Numeric_Kinds_Module
USE FileUnits_Asymptotic_Module
USE matrixlabels_Module      ! Use matrixlabels module for printing quantum labels
USE TotalEng_Module
USE Numbers_Module
USE Convrsns_Module
IMPLICIT NONE
CHARACTER(LEN=21), PARAMETER:: ProcName='Eigenphase'
CHARACTER(LEN=6) Print_Flag
!
! Author: Gregory A. Parker, Department of Physics and Astronomy, University of Oklahoma
!
! Calculates Eigenphases
!
! Required Input <=====
!     NOpen         Size of ALL matrices (NOpen by NOpen)
!     EDeriv=True   If ernergy derivatives of the K-Matrix exist
!     K_Matx        K-Matrix
!     KE_Matx       Energy derivative of the K-Matrix
!
! On RETURN  =====>

!     Phase         Eigenphases
!     Phase_Tran    Transformation to diagonalize the K-Matrix and the S-Matrix
!     PhaseE        Energy derivative of the eigenphases
!     PhaseE_Tran   Transformation to diagonalize the energy derivative of the K-Matrix.
!
! This routine is called by:
!      Asymptotic
! This routine calls:
!      DSYEV   a LAPACK Routine
!      Matrix_Out
!
INTEGER, INTENT(IN):: NOpen
LOGICAL, INTENT(IN):: EDeriv
INTEGER IDiag, LWork, Info
REAL(Kind=WP_Kind), INTENT(IN)::  K_Matx(NOpen,NOpen), KE_Matx(NOpen,NOpen)
REAL(Kind=WP_Kind), INTENT(OUT):: Phase(NOpen),  Phase_Tran(NOpen,NOpen)
REAL(Kind=WP_Kind), INTENT(OUT):: PhaseE(NOpen), PhaseE_Tran(NOpen,NOpen)
REAL(Kind=WP_Kind) K, KSQ
REAL(Kind=WP_Kind), ALLOCATABLE:: SRealDiag(:), SImagDiag(:), Work(:), Phase_Save(:)

WRITE(Out_Unit,*)
WRITE(Out_Unit,*)
WRITE(Out_Unit,*)
WRITE(Out_Unit,*)'Entering:', ProcName
CALL PoptAsy(ProcName, Print_Flag)

LWork=MAX(1,3*NOpen-1)
ALLOCATE(Work(LWork),Phase_Save(NOpen))
ALLOCATE(SRealDiag(NOpen), SImagDiag(NOpen))


Phase_Tran=K_Matx
CALL DSYEV('V','U', NOpen, Phase_Tran, NOpen, Phase, Work, LWork, Info)                ! Diagonalize K-Matrix
LabelRows=.False.
LabelColumns=.False.
CALL Matrix_Out(Phase, 1, Nopen, 'Phase', 'Eigenvalues of the K-Matrix', Print_Flag)   ! Print Eigenvalues of K-Matrix

!
! Calculate Eigenphases.
Phase_Save=Phase
DO IDiag=1,NOpen
   K=Phase(IDiag)
   KSQ=K*K
   SRealDiag(IDiag)=(1.d0-KSQ)/(1.d0+KSQ)
   SImagDiag(IDiag)=2.d0*K/(1.d0+KSQ)
   Phase(IDiag)=ATAN2(SImagDiag(IDiag),SRealDiag(IDiag))*0.5d0
ENDDO

!
! Store Eigenphases for Plotting
WRITE(eigunit,*)' $ready    Eigenphases'
WRITE(eigunit,*)NOpen
WRITE(eigunit,'(e14.7)')(Phase(IDiag),IDiag=1,NOpen)
WRITE(eigunit,*)' $CONTINUE'
IF(IthEnergy==1)WRITE(PhaseEPltUnit,'(A10,2A12,901(I12,","))')"IthEnergy","Etot(Ha)","Etot(ev)",(Idiag,idiag=1,nopen)
WRITE(PhasePltUnit,'(1x,I5,",",901(es12.5,","))')IthEnergy,Etot,Etot*autoev,(Phase(IDiag),IDiag=1,NOpen)

! Print Eigenphases
CALL Matrix_Out(Phase, 1, NOpen, 'Phase', 'Eigenphases', Print_Flag)
LabelColumns=.True.
! Print Transformation
CALL Matrix_Out(Phase_Tran, NOpen, NOpen, 'Phase_Tran', 'Transformation to diagonalize the S-Matrix',  Print_Flag)
Phase_Tran=Phase_Tran*Phase_Tran
! Print Phase Probabilities
CALL Matrix_Out(Phase_Tran, NOpen, NOpen, 'Phase_Prob', 'Eigenphase transformation probabilities',  Print_Flag)

IF(EDeriv)THEN
   PhaseE_Tran=KE_Matx
   ! Diagonalize energy derivative of K-Matrix
   CALL DSYEV('V','U', NOpen, PhaseE_Tran, NOpen, PhaseE, Work, LWork, Info)
   LabelColumns=.False.
   ! Print Eigenvalues of energy derivative
   CALL Matrix_Out(PhaseE, 1, NOpen, 'PhaseE', 'Eigenvalues of d(K-Matrix)/dE', Print_Flag)
   ! Calculate Eigenphases derivatives
   DO IDiag=1,NOpen
      K=Phase_Save(IDiag)
      KSQ=K*K
      PhaseE(IDiag)=PhaseE(IDiag)/(1.d0+KSQ)
   ENDDO
   ! Store Eigenphase derivatives for Plotting
   WRITE(derunit,*)' $ready    Energy Derivative of Eigenphases'
   WRITE(derunit,*)NOpen
   WRITE(derunit,'(e14.7)')(PhaseE(IDiag),IDiag=1,NOpen)
   WRITE(derunit,*)'0$CONTINUE'
   IF(IthEnergy==1)WRITE(PhaseEPltUnit,'(A10,2A12,901(I12,","))')"IthEnergy","Etot(Ha)","Etot(ev)",(Idiag,idiag=1,nopen)
   WRITE(PhaseEPltUnit,'(1x,I5,",",901(es12.5,","))')IthEnergy,Etot,Etot*autoev,(PhaseE(IDiag),IDiag=1,NOpen)
   ! Print eigenphase derivative
   CALL Matrix_Out(PhaseE, 1, NOpen, 'PhaseE', 'Energy derivative of Eigenphases', Print_Flag)
   LabelColumns=.True.
   ! Print Transformation
   CALL Matrix_Out(PhaseE_Tran, NOpen, NOpen, 'PhaseE_Tran', 'Transformation to diagonalize d(S-Matrix)/dE', Print_Flag)
   PhaseE_Tran=PhaseE_Tran*PhaseE_Tran
   CALL Matrix_Out(PhaseE_Tran, NOpen, NOpen, 'PhaseE_Prob', 'd(Eigenphase)/dE transformation probabilities', Print_Flag)
   ! Print PhaseE Probabilities
ENDIF

LabelRows=.True.
WRITE(Out_Unit,*)'Deallocate Temporary Storage in eigenphase'
DEALLOCATE(SRealDiag, SImagDiag, Work, Phase_Save)

WRITE(Out_Unit,*)'Leaving:', ProcName
RETURN
ENDSUBROUTINE Eigenphase
