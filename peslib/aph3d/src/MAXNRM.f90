SUBROUTINE MAXNRM(A, X, D, N, NORM)
USE Numeric_Kinds_Module
!
!   $RCSfile: maxnrm.f,v $ $Revision: 1.3 $
!   $Date: 89/07/28 09:57:22 $
!   $State: Stable $
!
!
! **********************************************************************
!
!     PURPOSE - (VER = 1)
!
!         COMPUTES THE MAXIMUM NORM OF A MATRIX STORED WITH THE
!         HELP OF THE D POINTER VECTOR.
!
!
!     INPUT PARAMETERS -
!
!         A = THE MATRIX.
!         D = POINTER VECTOR TO DIAGONAL ELEMENTS IN A.
!         N = DIMENSION OF A MATRIX.
!         X = A SCRATCH ARRAY OF LENGTH N.
!
!
!     OUTPUT PARAMETERS -
!
!         NORM = MAX ABS(A(I,1))+, ..., +ABS(A(I,N)).
!                 I
!
! **********************************************************************
!
REAL(Kind=WP_Kind)  A(1), DOT, NORM, X(1)
INTEGER   D(1), DI, DIM1, I, J, J2, LENM1, N, ROWNO, TOP
!
!     *******************************************
!     FOR I = 1 TO N DO
!       X(I) = ABS(A(I, 1)) + ... + ABS(A(I, N)).
!     *******************************************
DIM1 = 0
DO 30 I = 1, N
  TOP = DIM1 + 1
  DI = D(I)
  LENM1 = DI - TOP
!
  DOT = 0.0D0
  DO 10 J = TOP, DI
    DOT = DOT + ABS(A(J))
10        CONTINUE
  X(I) = DOT
!
  DIM1 = DI
  IF(.NOT. LENM1 > 0) GOTO 30
    ROWNO = I - LENM1
    J2 = I - 1
    DO 20 J = ROWNO, J2
      X(J) = X(J) + ABS(A(TOP))
      TOP = TOP + 1
20          CONTINUE
!
30      CONTINUE
!
!     *********************************
!     COMPUTE MAX(X(I)), I = 1, ..., N.
!     *********************************
NORM = X(1)
DO 40 I = 1, N
  NORM = MAX(NORM, X(I))
40      CONTINUE
!
RETURN
END
