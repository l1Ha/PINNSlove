LOGICAL FUNCTION IMTQL2(ALPHA, BETA, D, E, Z, N, NM, FINAL)
USE Numeric_Kinds_Module
!
!   $RCSfile: imtql2.f,v $ $Revision: 1.3 $
!   $Date: 89/07/28 09:57:13 $
!   $State: Stable $
!
!
! **********************************************************************
!
!     PURPOSE -
!
!         A SOMEWHAT MODIFIED VERSION OF THE EISPACK SUBROUTINE
!         IMTQL2. SEE B.S. GARBOW AND J.J. DONGARRA, PATH CHART AND
!         DOCUMENTATION FOR THE EISPACK PACKAGE OF MATRIX EIGENSYSTEM
!         ROUTINES, ARGONNE NATIONAL LABORATORY, 1975.
!
!         THE MAJOR MODIFICATIONS ARE,
!
!                 DIFFERENT NUMBER OF DUMMY PARAMETERS.
!                 INITIALIZATION OF D, E, AND, Z INSIDE THE PROCEDURE.
!                 POSSIBILTY TO COMPUTE EIGENVALUES AND ONLY TOP AND
!                 BOTTOM ELEMENTS OF Z (THE EIGENVECTORS).
!                 NO SORTING IS MADE.
!                 SINGLE PRECISION.
!
!         THE LOOP WHERE THE COMPUTATIONS ARE MADE (DO 240 ...)
!         IS HOWEVER NEARLY INTACT. IN THIS CASE WE WILL DESCRIBE
!         THE PARAMETERS SINCE THE NAMES ARE DIFFERENT FROM WHAT WE
!         ARE USED TO.
!
!         IMTQL2 IS USED AS FOLLOWS IN THE REFERENCING ROUTINE,
!
!
!         IF(.NOT. IMTQL2(W(ALPHA), W(BETA), W(NU), W(SCR), W(S),
!        +                P, PMAX, FINAL)) GOTO 8888
!
!
!     INPUT PARAMETERS -
!
!         ALPHA = DIAGONAL IN THE TRIDIAGONAL MATRIX T.
!         BETA  = SUBDIAGONAL IN T.
!         E     = IS A WORKING AREA.
!         N     = THE DIMENSION OF T.
!         NM    = ROW DIMENSION OF Z, THE EIGENVECTOR MATRIX.
!         FINAL = IF TRUE, COMPUTE THE WHOLE EIGENVECTORS,
!                 IF FALSE, COMPUTE ONLY THE TOP AND BOTTOM ELEMENTS.
!
!
!     OUTPUT PARAMETERS -
!
!         D = CONTAINS THE EIGENVALUES.
!         Z = CONTAINS THE EIGENVECTORS OR THE END ELEMENTS.
!
!
!     NOTE. IN THE ORIGINAL VERSION, IMTQL2 OVERWRITES ALPHA WITH
!     THE EIGENVALUES AND DESTROYS BETA. Z SHOULD FURTHERMORE BE SET TO
!     THE IDENTITY MATRIX. WE THEREFORE COPY ALPHA TO D, BETA TO E,
!     AND Z IS SET TO I (OR THE FIRST AND LAST ROWS OF I) BEFORE
!     ANY COMPUTATION TAKES PLACE.
!
! **********************************************************************
!
!
INTEGER   NM
REAL(Kind=WP_Kind) ALPHA(1), B, BETA(1), C, D(1), E(1), F, FACTOR, G,  SRELPR, P, R, RDUMP, S, SECOND, ST, TIME, TOLBPI, &
               TOLLDL, TOLPDM, TOLS1I, TOLZBT, TOLZNU, Z(NM, 1)
INTEGER   COUNT, DUMMY, ERRNO, I, IDUMP, II, J, K, L, M, MML, N,  NBADMU, NMXRST, NOR, STEP
LOGICAL   FA, FINAL, T
COMMON   /STLMER/ RDUMP, ERRNO, IDUMP(2)
COMMON   /STLMST/ TIME(24), COUNT(24), NBADMU, NMXRST, DUMMY
COMMON   /STLMTF/ T, FA
COMMON   /STLMTL/ SRELPR, TOLBPI, TOLS1I, FACTOR, TOLPDM, TOLZBT, TOLZNU, TOLLDL(3)
!
DATA      NOR  / 6 /
!
ST = SECOND()
COUNT(NOR) = COUNT(NOR) + 1
IMTQL2 = T
!
!     *****
!     COPY.
!     *****
DO 10 I = 1, N
  D(I) = ALPHA(I)
  E(I) = BETA(I)
10      CONTINUE
E(N) = 0.0D0
!
!     *************
!     Z = IDENTITY.
!     *************
IF(.NOT. FINAL) GOTO 40
  STEP = 1
  DO 30 I = 1, N
    DO 20 J = 1, N
      Z(I, J) = 0.0D0
20          CONTINUE
    Z(I, I) = 1.0D0
30        CONTINUE
  GOTO 60
!
!     **********************************
!     Z = IDENTITY, FIRST AND LAST ROWS.
!     **********************************
40    STEP = N - 1
DO 50 I = 1, N
  Z(1, I) = 0.0D0
  Z(N, I) = 0.0D0
50      CONTINUE
Z(1, 1) = 1.0D0
Z(N, N) = 1.0D0
!
60    IF(N == 1) GOTO 9999
!
DO 240 L = 1, N
   J = 0
!        ***************************************
!        LOOK FOR SMALL SUB - DIAGONAL ELEMENTS.
!        ***************************************
105      DO 110 M = L, N
      IF(M == N) GOTO 120
      IF(ABS(E(M)) <= SRELPR * (ABS(D(M)) + ABS(D(M + 1))))GOTO 120
110      CONTINUE
!
120      P = D(L)
   IF(M == L) GOTO 240
   IF(.NOT. J == 30) GOTO 121
!        ***************
!        ERROR HANDLING.
!        ***************
    IDUMP(1) = L
    CALL ERROR(NOR, 1)
    GOTO 8888
!
121      J = J + 1
!        ***********
!        FORM SHIFT.
!        ***********
   G = (D(L + 1) - P) / (2.0D0 * E(L))
   R = SQRT(G * G + 1.0D0)
   G = D(M) - P + E(L) / (G + SIGN(R, G))
   S = 1.0D0
   C = 1.0D0
   P = 0.0D0
   MML = M - L
!        **********************************
!        FOR I = M - 1 STEP  - 1 DO  -  -
!        **********************************
   DO 200 II = 1, MML
      I = M - II
      F = S * E(I)
      B = C * E(I)
      IF(ABS(F) < ABS(G)) GOTO 150
      C = G / F
      R = SQRT(C * C + 1.0D0)
      E(I + 1) = F * R
      S = 1.0D0 / R
      C = C * S
      GOTO 160
150         S = F / G
      R = SQRT(S * S + 1.0D0)
      E(I + 1) = G * R
      C = 1.0D0 / R
      S = S * C
160         G = D(I + 1) - P
      R = (D(I) - G) * S + 2.0D0 * C * B
      P = S * R
      D(I + 1) = G + P
      G = C * R - B
!           *****************************************************
!           FORM VECTOR.
!           STEP = 1, COMPUTE THE WHOLE VECTORS.
!           STEP = N-1, COMPUTE ONLY THE TOP AND BOTTOM ELEMENTS.
!           *****************************************************
      DO 180 K = 1, N, STEP
         F = Z(K, I + 1)
         Z(K, I + 1) = S * Z(K, I) + C * F
         Z(K, I) = C * Z(K, I) - S * F
180         CONTINUE
!
200      CONTINUE
!
   D(L) = D(L) - P
   E(L) = G
   E(M) = 0.0D0
   GOTO 105
240   CONTINUE
!
GOTO 9999
8888  IMTQL2 = FA
CALL ERROR(NOR, NOR)
!
9999  TIME(NOR) = TIME(NOR) + SECOND() - ST
RETURN
END
