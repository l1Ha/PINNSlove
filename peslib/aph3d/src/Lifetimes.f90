SUBROUTINE Lifetimes(NOpen, S_Real, S_Imag, SE_Real, SE_Imag,  Q_Real, Q_Imag)
USE Numeric_Kinds_Module       ! Set numeric kinds
USE FileUnits_Asymptotic_Module           ! Fileunit numbers
USE Matrixlabels_Module        ! Use matrixlabels module for printing quantum labels
USE Physical_Constants_Module  ! Contains physical constants
USE Convrsns_Module
USE TotalEng_Module
IMPLICIT NONE
CHARACTER(LEN=9), PARAMETER:: ProcName='Lifetimes'
CHARACTER(LEN=6) Print_Flag
!
!     Author: Gregory A. Parker, Department of Physics and Astronomy, University of Oklahoma
!
!     This routine determines the lifetime matrix Q as:
!     Q = -i*hbar*(dS/dE)*S(transpose conjugate)
!          where S is the S-Matrix and the derivative is with
!          respect to the total energy (E)
!
!     References
!     Z. Darakjian, E.F. Hayes, G.A. Parker, E.A. Butcher and J.D. Kress JCP 95, 2516-2522 (1991)
!     F.T. Smith, Phys. Rev. 118, 349 (1960) [Original definition of Lifetime Matrix]
!
!     Required Input <=====
!     NOpen           Size of ALL matrices (NOpen by NOpen)
!     S_Real      Real part of the S-Matrix
!     S_Imag      Imaginary part of the S-Matrix
!     SE_Real     Energy derivative of the real part of the S-Matrix
!     SE_Imag     Energy derivative of the imaginary part of the S-Matrix
!                 The units of Energy must be Hartree Atomic Units.
!
!     On return  =====>
!     Q_Real      Real part of the Lifetime Matrix (Q)
!     Q_Imag      Imaginary part of the Lifetime Matrix (Q)
!                 Units of the lifetime matrix are Femto-seconds
!
!     This routine is called by:
!          Asymptotic
!     This routine calls:
!          Matrix_Out
!
INTEGER, INTENT(IN)                :: NOpen
INTEGER                            :: Info, LWork, I, J
REAL(KIND=WP_Kind), INTENT(IN)     :: S_Real(NOpen,NOpen),  S_Imag(NOpen,NOpen)
REAL(KIND=WP_Kind), INTENT(IN)     :: SE_Real(NOpen,NOpen), SE_Imag(NOpen,NOpen)
REAL(KIND=WP_Kind), INTENT(OUT)    :: Q_Real(NOpen,NOpen),  Q_Imag(NOpen,NOpen)
REAL(KIND=WP_Kind), ALLOCATABLE    :: RWork(:,:), Eigenlife(:)
COMPLEX(KIND=WP_Kind), ALLOCATABLE :: QTran(:,:), Work(:,:)


WRITE(Out_Unit,*)
WRITE(Out_Unit,*)
WRITE(Out_Unit,*)
WRITE(Out_Unit,*)'Entering:', ProcName
CALL PoptAsy(ProcName, Print_Flag)

LWORK=NOpen*NOpen
ALLOCATE(QTran(NOpen,NOpen),Work(NOpen,NOpen),RWork(NOpen,NOpen),Eigenlife(NOpen))
!
!     real part of Q = SE_Imag*S_Real(trans) - SE_Real*S_Imag(trans)
Q_Real= MATMUL(SE_Imag,TRANSPOSE(S_Real))-MATMUL(SE_Real,TRANSPOSE(S_Imag))    ! Calculate real part of Lifetime matrix

!
!     imag part of Q = -[SE_Real*S_Real(trans) + SE_Imag*S_Imag(trans)]
Q_Imag=-MATMUL(SE_Real,TRANSPOSE(S_Real))-MATMUL(SE_Imag,TRANSPOSE(S_Imag))    ! Calculate real part of Lifetime matrix

!
!    Convert to desired units (In this case femtoseconds)
WRITE(Out_Unit,*)'Time_Conv=',Time_Conv
Q_Real = mhbar*Q_Real
Q_Imag = mhbar*Q_Imag

!
!     Print out part of the Liftime Matrix (Q)
CALL Matrix_Out (Q_Real, NOpen, NOpen,'Q_Real','Real Part of Lifetime Matrix',  Print_Flag)         ! Print real part
CALL Matrix_Out (Q_Imag, NOpen, NOpen,'Q_Imag','Imaginary Part of Lifetime Matrix',   Print_Flag)   ! Print imaginary part

!
! Find the unitarity transformations (QTran) which diagonalizes the Q-Matrix.
!      The eigenvalues are the eigenlifetimes.
! QProb=square magnitude of (QTran) which gives the probability that a particular
!      quantum state contributes to a given Eigenlifetime
QTran=CMPLX(Q_Real,Q_Imag,KIND=WP_Kind)
CALL ZHEEV('V', 'U', Nopen, QTran, NOpen, Eigenlife, Work, LWork, RWork, Info)
IF(Info/=0)THEN
   WRITE(OUT_Unit,*)'Error in Lifetimes'
   !STOP 'Lifetimes'  GregParker TempMod
ENDIF

LabelRows=.False.
LabelColumns=.False.
CALL Matrix_Out(Eigenlife, 1, NOpen, 'Eigenlife', 'Eigen-Lifetimes', Print_Flag)                   ! Print eigenlifetimes
LabelColumns=.True.
IF(IthEnergy==1)WRITE(PhaseEPltUnit,'(A10,2A12,901(I12,","))')"IthEnergy","Etot(Ha)","Etot(ev)",(I,I=1,nopen)
WRITE(EigLifePltUnit,'(1x,I5,",",901(es12.5,","))')IthEnergy,Etot,Etot*autoev,(Eigenlife(I),I=1, NOpen)

!
! Store the real part of the transformation in a real array for printing
DO I=1,NOpen
   DO J=1,NOpen
      RWork(I,J)=REAL(QTran(I,J))
   ENDDO
ENDDO

CALL Matrix_Out(RWork, NOpen, NOpen, 'QTran_Real', 'Real Part of QTran', Print_Flag)               ! Print Real(QTran)

!
! Store the imaginary part of the transformation in a real array for printing
DO I=1,NOpen
   DO J=1,NOpen
      RWork(I,J)=AIMAG(QTran(I,J))
   ENDDO
ENDDO

CALL Matrix_Out(RWork, NOpen, NOpen, 'QTran_Imag', 'Imaginary Part of QTran', Print_Flag)          ! Print Imag(QTran)

!
! Calculate magnitude squared of the transformation matrix.
! The elements are the probabilites that each state contributes to a particular eigenlifetime
DO I=1,NOpen
   DO J=1,NOpen
      RWork(I,J)=REAL(QTran(I,J))*REAL(QTran(I,J))+AIMAG(QTran(I,J))*AIMAG(QTran(I,J))
   ENDDO
ENDDO

CALL Matrix_Out(RWork, NOpen, NOpen, 'QProb', 'QTranProbability', Print_Flag)                      ! Print QTranProbability

LabelRows=.True.


DEALLOCATE(QTran, Work, RWork, Eigenlife)

WRITE(Out_Unit,*)'Leaving:', ProcName
RETURN
ENDSUBROUTINE Lifetimes
