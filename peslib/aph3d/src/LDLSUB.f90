SUBROUTINE LDLSUB(K, M, D, MXNORM, ONENRM, MU, N, RNEW, DIAGM, MEQI, LEN)
USE Numeric_Kinds_Module
!
!   $RCSfile: ldlsub.f,v $ $Revision: 1.3 $
!   $Date: 89/07/28 09:57:21 $
!   $State: Stable $
!
!
! **********************************************************************
!
!     PURPOSE - (VER = 1)
!
!         COMPUTES K-MU*M AND MAKES AN LDL(T)-DECOMPOSITION OF IT.
!         COMPUTES THE NUMBER OF NEGATIVE D(I,I)-ELEMENTS.
!
!
!     INPUT PARAMETERS -
!
!         D       THE POINTER VECTOR TO THE DIAGONAL ELEMENTS IN K.
!         DIAGM   PROFIL <= 2.
!         K       K MATRIX (STORED USING PROFILE STORAGE).
!         LEN     LOCAL ERROR NUMBER.
!         M       M MATRIX (STORED ACCORDING TO THE USER GUIDE).
!         MEQI    PROFIL == 3
!         N       DIMENSION OF PROBLEM.
!
!         MXNORM AND ONENRM ARE SCRATCH VECTORS OF LENGTH N.
!
!     OUTPUT PARAMETERS -
!
!         K       IS OVERWRITTEN AND CONTAINS L AND D, IN THE LDL(T)
!                 FACTORIZATION OF K-MU*M.
!
!     PLEASE SEE THE PROGRAMMERS GUIDE FOR INFORMATION ABOUT
!     PARAMETERS NOT EXPLAINED ABOVE, AND FOR MORE DETAILS ABOUT
!     THE FUNCTION OF THE ROUTINE.
!
!
! **********************************************************************
!
REAL(Kind=WP_Kind)  DIAGEL, FACTOR, K(1), LIJ, M(1), SRELPR, MU, MXNORM(1)
REAL(Kind=WP_Kind)  NORM, ONENRM(1), RDUMP, SCPR, SIJ, TOL, TOLBPI, TOLLDL
REAL(Kind=WP_Kind)  TOLPDM, TOLS1I, TOLZBT, TOLZNU
INTEGER D(1), DI, DIM1, DIMI, DJ, ERRNO, I, IDUMP, IJAD, IJAD1,  IJAD2
INTEGER J, LEN, LENM1, LOOPL, MXQIQJ, N, NM1, QI, RNEW
INTEGER STARTL, STARTS, TOP
LOGICAL   DIAGM, F, MEQI, T
COMMON   /STLMER/ RDUMP, ERRNO, IDUMP(2)
COMMON   /STLMTF/ T, F
COMMON   /STLMTL/ SRELPR, TOLBPI, TOLS1I, FACTOR, TOLPDM, TOLZBT,  TOLZNU, TOLLDL(3)
!
RNEW = 0
LEN = 0
!

!     ********************
!     FORM K = K - MU * M.
!     ********************
IF(.NOT. DIAGM) GOTO 40
  IF(.NOT. MEQI) GOTO 20
!         ******
!         M = I.
!         ******
    DO 10 I = 1, N
      DI = D(I)
      K(DI) = K(DI) - MU
10          CONTINUE
    GOTO 50
!
!       **********************************
!       M IS DIAGONAL, BUT NOT EQUAL TO I.
!       **********************************
20      DO 30 I = 1, N
    DI = D(I)
    K(DI) = K(DI) - MU * M(I)
30        CONTINUE
  GOTO 50
!
!     ****************************
!     M HAS THE SAME PROFILE AS K.
!     ****************************
40    CALL SUBV(K, M, MU, D(N))
!
!     *****************************
!     COMPUTE MAX NORM(K - MU * M),
!     AND THE FIRST TOLERANCE.
!     *****************************
50    CALL MAXNRM(K, MXNORM, D, N, NORM)
IF(NORM == 0.0D0) NORM = 1.0D0
TOL = TOLLDL(1) * NORM
!
DO 60 I = 1, N
  MXNORM(I) = 0.0D0
  ONENRM(I) = 0.0D0
60      CONTINUE
!
DIM1 = 0
!
!     *********************
!     FOR I = 1 TO N DO ...
!     *********************
DO 120 I = 1, N
  TOP = DIM1 + 1
  DI = D(I)
  LENM1 = DI - TOP
!       *****************************************
!       IF(QI <= I - 1)THEN ... ELSE GOTO 110
!       *****************************************
  IF(.NOT. LENM1 >= 1) GOTO 110
    QI = I - LENM1
    J = QI
    IJAD2 = DI - 1
!         ****************************************
!         IF(QI < I - 1)THEN ... ELSE GOTO 80
!         ****************************************
    IF(.NOT. LENM1 > 1) GOTO 80
      IJAD1 = TOP + 1
      DIMI = DI - I
!           ******************************
!           FOR J = QI + 1 TO I - 1 DO ...
!           ******************************
      DO 70 IJAD = IJAD1, IJAD2
        J = J + 1
        DJ = D(J)
!             **********************
!             MXQIQJ = MAX (QI, QJ).
!             **********************
        MXQIQJ = MAX0(QI, J - DJ + D(J - 1) + 1)
        LOOPL = J - MXQIQJ
        IF(.NOT. LOOPL >= 1) GOTO 70
!               *****************************
!               COMPUTE SUM L(J, U) * S(I, U),
!               FOR U = MAX(QI, QJ) TO J - 1.
!               *****************************
          STARTL = DJ - LOOPL
          STARTS = DIMI + MXQIQJ
          K(IJAD) = K(IJAD) - SCPR(K(STARTL), K(STARTS), LOOPL)
70            CONTINUE
!
    J = QI
80        DIAGEL = K(DI)
!         **************************
!         FOR J = QI TO I - 1 DO ...
!         **************************
    DO 100 IJAD = TOP, IJAD2
      DJ = D(J)
      SIJ = K(IJAD)
      LIJ = SIJ / K(DJ)
      ONENRM(J) = ONENRM(J) + ABS(LIJ)
      MXNORM(J) = MAX(MXNORM(J), ABS(SIJ))
      DIAGEL = DIAGEL - LIJ * SIJ
      K(IJAD) = LIJ
      J = J + 1
100         CONTINUE
!
    K(DI) = DIAGEL
!
110     DIAGEL = K(DI)
!       ************
!       UPDATE RNEW.
!       ************
  IF(DIAGEL < 0.0D0) RNEW = RNEW + 1
  DIM1 = DI
!
!       ***************************
!       IF D(J) IS TOO SMALL, STOP.
!       ***************************
  IF(.NOT. ABS(DIAGEL) <= TOL) GOTO 120
    LEN = -1
    IDUMP(1) = I
    IDUMP(2) = 1
    RDUMP = DIAGEL
    GOTO 9999
!
120     CONTINUE
!
IF(N == 1) GOTO 9999
!     ****************************
!     COMPUTE THE OTHER TOLERANCE.
!     ****************************
TOL = TOLLDL(2) * NORM
NM1 = N - 1
!     *************************
!     FOR I = 1 TO N - 1 DO ...
!     *************************
DO 130 I = 1, NM1
  IF(.NOT. ONENRM(I) * MXNORM(I) > TOL) GOTO 130
    LEN = -1
    IDUMP(1) = I
    IDUMP(2) = 1
    RDUMP = ONENRM(I) * MXNORM(I)
    GOTO 9999
!
130     CONTINUE
!
9999  RETURN
END
