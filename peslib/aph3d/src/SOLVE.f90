SUBROUTINE SOLVE(X, B, K, D, N)
USE Numeric_Kinds_Module
!
!   $RCSfile: solve.f,v $ $Revision: 1.3 $
!   $Date: 89/07/28 09:57:42 $
!   $State: Stable $
!
!
! **********************************************************************
!
!     PURPOSE - (VER = 1)
!
!         SOLVES (LDL(T))*X=B WITH THE USUAL METHOD, VIZ.
!
!         1) SOLVE L*Z=B      (FORWARD SUBSTITUTION)
!         2) COMPUTE U=(D**(-1))*Z
!         3) SOLVE L(T)*X=U   (BACK SUBSTITUTION)
!
!         IN THIS ROUTINE Z,U, AND X OCCUPY THE SAME MEMORY.
!
!
!     INPUT PARAMETERS -
!
!         B = KNOWN VECTOR.
!         D = POINTER VECTOR TO THE DIAGONAL ELEMENTS IN L.
!         K = CONTAINS THE STRICTLY LOWER TRIANGLE OF L, AND
!             THE DIAGONAL, D, OF LDL(T) (STORED IN THE DIAGONAL OF K).
!         N = DIMENSION OF L MATRIX.
!
!
!     OUTPUT PARAMETERS -
!
!         X = SOLUTION.
!
! **********************************************************************
!
REAL(Kind=WP_Kind)  B(1), K(1), SCPR, X(1), XNDECR
INTEGER   D(1), DI, DIM1, I, LENM1, N, NDECR, ROWNO, TOP
!
!     *********************
!     FORWARD SUBSTITUTION.
!     *********************
DIM1 = 0
DO 20 I = 1, N
  DI = D(I)
  TOP = DIM1 + 1
  LENM1 = DI - TOP
  DIM1 = DI
  ROWNO = I - LENM1
  IF(.NOT. LENM1 > 0) GOTO 10
    X(I) = B(I) - SCPR(K(TOP), X(ROWNO), LENM1)
    GOTO 20
10      X(I) = B(I)
20      CONTINUE
!
!     **********
!     D-INVERSE.
!     **********
DO 40 I = 1, N
  DI = D(I)
  X(I) = X(I) / K(DI)
40      CONTINUE
!
IF(N == 1) GOTO 9999
!
!     ******************
!     BACK SUBSTITUTION.
!     ******************
NDECR = N
DO 50 I = 2, N
  TOP = D(NDECR - 1) + 1
  LENM1 = D(NDECR) - TOP
  ROWNO = NDECR - LENM1
  XNDECR = X(NDECR)
  IF(LENM1 > 0) CALL SUBV(X(ROWNO), K(TOP), XNDECR, LENM1)
  NDECR = NDECR - 1
50      CONTINUE
!
9999  RETURN
END
