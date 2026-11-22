SUBROUTINE Numerov_Matrix(N, hstep, usys2, X_Values, V, Eig_Values, Eig_Vectors, amat, bmat, neig, MaxVibFull)
USE Numeric_Kinds_Module
USE Numbers_Module
USE CriticalDist_Module
USE FileUnits_Module
USE FileUnits_OneDim_Module
USE QstateDAF_Module
USE OneDim_Parms_Module
USE NArran_Module
USE DiatomicPot_Module
USE Convrsns_Module
!Matrix Numerov
!See Am. J. Phys. 80 1017-1019 (2012) : Mohandas Pillai, Joshua Goglio and Thad Walker
USE Numeric_Kinds_Module
IMPLICIT NONE
CHARACTER(LEN=1), PARAMETER:: Jobz='V', UpLo='L'
INTEGER, PARAMETER:: lwork=100000
INTEGER N, i, j, info, NBound, neig, MaxVibFull
! PORT FIX (2026-10-01): S_work must be double precision; it was INTEGER, so
! LAPACK wrote 2*lwork doubles into an lwork-integer buffer and corrupted
! memory (SIGSEGV inside DSYSV/DSYCONV).
INTEGER IPIV(lwork)
REAL(KIND=WP_Kind) S_work(lwork)
REAL(KIND=WP_Kind) hstep, usys2
REAL(KIND=WP_Kind) Amat(N,N), Bmat(N,N), Eig_Values(N), Eig_Vectors(N,N), X_Values(N), V(N)
!REAL(KIND=WP_Kind), PARAMETER::  one=1.0_WP_Kind, two=2.0_WP_Kind, ten=10.0_WP_Kind, twelve=12.0_WP_Kind

!create A and B Matricex See equations 3-5 in reference
Amat=0.d0
Bmat=0.d0
DO I=1,N
    Amat(I,I)=-two
    Bmat(I,I)=ten/twelve
    IF(I>1)THEN
        Amat(I-1,I)=one
        Bmat(I-1,I)=one/twelve
    ENDIF
    IF(I<N)THEN
        Amat(I+1,I)=one
        Bmat(I+1,I)=one/twelve
    ENDIF
ENDDO

!Create the Kinetic Energy Matrix
Amat=-Amat/(usys2*hstep**2)
CALL DSYSV('U', N, N, Bmat, N, IPIV, Amat, N, s_work, lwork, INFO)

!Calculate the Hamiltonian Matrix H=T+V
DO i=1,N
    Amat(i,i)=Amat(i,i)+V(i)
ENDDO

!Find the Eigenvalues and Eigenvectors (Eigenfunctions evaluated at the Grid Points) of the Hamiltonian
CALL dsyev(Jobz, Uplo, N, Amat, N, Eig_Values, s_work, lwork, info)

!Write the eigenvalues
OPEN(Unit=151,File=TRIM(DiatomicPot)//'_'//"Eigenvalues_Numerov.csv")
WRITE(151,*)TRIM(DiatomicPot)
WRITE(151,*)"Info=",Info
WRITE(151,*)"Numerov Eigenvalues"
WRITE(131,*)"State,  Eig(Hartree),  Eig(cm^-1), Type"
DO i=1,N
   IF(V(N)>=Eig_Values(i))THEN
       WRITE(151,'(I5,",",ES22.15,",",ES22.15,",",A)')i,Eig_Values(i),Eig_Values(i)*1822.0,"  Bound"
       NBound=i
   ELSE
       WRITE(151,'(I5,",",ES22.15,",",ES22.15,",",A)')i,Eig_Values(i),Eig_Values(i)*1822.0,"  Continuum"
   ENDIF
ENDDO
CLOSE(Unit=151)

!Write the eigenvectors
OPEN(Unit=151,File=TRIM(DiatomicPot)//'_'//"Eigenfunctions_Numerov.csv")
WRITE(151,*)TRIM(DiatomicPot)
WRITE(151,*)"Info=",Info
WRITE(151,*)"Numerov Eigenvectors"
DO i=1,N
   WRITE(151,'(201(ES14.5,","))')X_Values(i),(Amat(i,j),j=1,NBound)
ENDDO
CLOSE(Unit=151)

RETURN
ENDSUBROUTINE Numerov_Matrix