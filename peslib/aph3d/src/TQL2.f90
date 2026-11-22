SUBROUTINE TQL2(NM,N,D,E,Z,IERR)
USE Numeric_Kinds_Module
USE Numbers_Module
!   BEGIN PROLOGUE  TQL2
!   DATE WRITTEN   760101   (YYMMDD)
!   REVISION DATE  861211   (YYMMDD)
!   CATEGORY NO.  D4A5,D4C2A
!   KEYWORDS  LIBRARY=SLATEC(EISPACK),TYPE=SINGLE PRECISION(TQL2-S),
!             EIGENVALUES,EIGENVECTORS
!   AUTHOR  SMITH, B. T., ET AL.
!   PURPOSE  Compute eigenvalues and eigenvectors of symmetric
!            tridiagonal matrix.
!   DESCRIPTION
!
!     This SUBROUTINE is a translation of the ALGOL procedure TQL2,
!     NUM. MATH. 11, 293-306(1968) by Bowdler, Martin, Reinsch, and
!     Wilkinson.
!     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 227-240(1971).
!
!     This SUBROUTINE finds the eigenvalues and eigenvectors
!     of a SYMMETRIC TRIDIAGONAL matrix by the QL method.
!     The eigenvectors of a FULL SYMMETRIC matrix can also
!     be found IF  TRED2  has been used to reduce this
!     full matrix to tridiagonal form.
!
!     On Input
!
!        NM must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        N is the order of the matrix.
!
!        D contains the diagonal elements of the input matrix.
!
!        E contains the subdiagonal elements of the input matrix
!          in its last N-1 positions.  E(1) is arbitrary.
!
!        Z contains the transformation matrix produced in the
!          reduction by  TRED2, IF performed.  IF the eigenvectors
!          of the tridiagonal matrix are desired, Z must contain
!          the identity matrix.
!
!      On Output
!
!        D contains the eigenvalues in ascending order.  IF an
!          error exit is made, the eigenvalues are correct but
!          unordered for indices 1,2,...,IERR-1.
!
!        E has been destroyed.
!
!        Z contains orthonormal eigenvectors of the symmetric
!          tridiagonal (or full) matrix.  IF an error exit is made,
!          Z contains the eigenvectors associated with the stored
!          eigenvalues.
!
!        IERR is set to
!          Zero       for normal RETURN,
!          J          IF the J-th eigenvalue has not been
!                     determined after 30 iterations.
!
!     Calls PYTHAG(A,B) for sqrt(A**2 + B**2).
!
!     Questions and comments should be directed to B. S. Garbow,
!     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
!     ------------------------------------------------------------------
!   REFERENCES  B. T. SMITH, J. M. BOYLE, J. J. DONGARRA, B. S. GARBOW,
!                 Y. IKEBE, V. C. KLEMA, C. B. MOLER, *MATRIX EIGEN-
!                 SYSTEM ROUTINES - EISPACK GUIDE*, SPRINGER-VERLAG,
!                 1976.
!   ROUTINES CALLED  PYTHAG
!   END PROLOGUE  TQL2
!
INTEGER I,J,K,L,M,N,II,L1,L2,NM,MML,IERR
REAL(Kind=WP_Kind) D(N),E(N),Z(NM,N)
REAL(Kind=WP_Kind) B,C,C2,C3,DL1,EL1,F,G,H,P,R,S,S2
REAL(Kind=WP_Kind) PYTHAG
!
!   FIRST EXECUTABLE STATEMENT  TQL2
IF(n==0)RETURN   !TMPMODGregParker
IERR = 0
IF(N == 1) GO TO 1001
!
DO I = 2, N
   E(I-1) = E(I)
ENDDO
!
F = 0.0D0
B = 0.0D0
E(N) = 0.0D0
!
DO 240 L = 1, N
   J = 0
   H = ABS(D(L)) + ABS(E(L))
   IF(B < H) B = H
!     .......... LOOK FOR SMALL SUB-DIAGONAL ELEMENT ..........
   DO 110 M = L, N
      IF(B + ABS(E(M)) == B) GO TO 120
!     .......... E(N) IS ALWAYS ZERO, SO THERE IS NO EXIT
!                THROUGH THE BOTTOM OF THE LOOP ..........
  110    CONTINUE
!
  120    IF(M == L) GO TO 220
  130    IF(J == 30) GO TO 1000
   J = J + 1
!     .......... FORM SHIFT ..........
   L1 = L + 1
   L2 = L1 + 1
   G = D(L)
   P = (D(L1) - G) / (Two * E(L))
   R = PYTHAG(P,One)
   D(L) = E(L) / (P + SIGN(R,P))
   D(L1) = E(L) * (P + SIGN(R,P))
   DL1 = D(L1)
   H = G - D(L)
   IF(L2 > N) GO TO 145
!
   DO I = L2, N
      D(I) = D(I) - H
   ENDDO
!
  145    F = F + H
!     .......... QL TRANSFORMATION ..........
   P = D(M)
   C = 1.0D0
   C2 = C
   EL1 = E(L1)
   S = 0.0D0
   MML = M - L
!     .......... FOR I=M-1 STEP -1 UNTIL L DO -- ..........
   DO 200 II = 1, MML
      C3 = C2
      C2 = C
      S2 = S
      I = M - II
      G = C * E(I)
      H = C * P
      IF(ABS(P) < ABS(E(I))) GO TO 150
      C = E(I) / P
      R = SQRT(C*C+1.0D0)
      E(I+1) = S * P * R
      S = C / R
      C = 1.0D0 / R
      GO TO 160
  150       C = P / E(I)
      R = SQRT(C*C+1.0D0)
      E(I+1) = S * E(I) * R
      S = 1.0D0 / R
      C = C * S
  160       P = C * D(I) - S * G
      D(I+1) = H + S * (C * G + S * D(I))
!     .......... FORM VECTOR ..........
      DO 180 K = 1, N
         H = Z(K,I+1)
         Z(K,I+1) = S * Z(K,I) + C * H
         Z(K,I) = C * Z(K,I) - S * H
  180       CONTINUE
!
  200    CONTINUE
!
   P = -S * S2 * C3 * EL1 * E(L) / DL1
   E(L) = S * P
   D(L) = C * P
   IF(B + ABS(E(L)) > B) GO TO 130
  220    D(L) = D(L) + F
  240 CONTINUE
!     .......... ORDER EIGENVALUES AND EIGENVECTORS ..........
DO 300 II = 2, N
   I = II - 1
   K = I
   P = D(I)
!
   DO 260 J = II, N
      IF(D(J) >= P) GO TO 260
      K = J
      P = D(J)
  260    CONTINUE
!
   IF(K == I) GO TO 300
   D(K) = D(I)
   D(I) = P
!
   DO 280 J = 1, N
      P = Z(J,I)
      Z(J,I) = Z(J,K)
      Z(J,K) = P
  280    CONTINUE
!
  300 CONTINUE
!
GO TO 1001
!     .......... SET ERROR -- NO CONVERGENCE TO AN
!                EIGENVALUE AFTER 30 ITERATIONS ..........
 1000 IERR = L
 1001 RETURN
END
