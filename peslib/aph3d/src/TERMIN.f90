LOGICAL FUNCTION TERMIN(CAUSE, W, IW, BG, TCONV, NLEFT, ERRNO)
USE Numeric_Kinds_Module
!
!   $RCSfile: termin.f,v $ $Revision: 1.3 $
!   $Date: 89/07/28 09:57:48 $
!   $State: Stable $
!
!
! **********************************************************************
!
!     PURPOSE -
!
!         THIS ROUTINE CHECKS IF WE HAVE COMPUTED THE
!         REQUESTED EIGENPAIRS OR NOT.
!
!     INPUT PARAMETERS -
!
!         CAUSE = 1, DECIDE IF WE SHOULD CONTINUE OR NOT.
!                 2, WE HAVE REACHED THE LAST EIGENVALUE.
!                 3, AN ERROR HAS OCCURRED.
!
!     OUTPUT PARAMETERS -
!
!         TERMIN  = FALSE, IF WE SHOULD TAKE ANOTHER SHIFT,
!                   AND TRUE OTHERWISE.
!
!     PLEASE SEE THE PROGRAMMERS AND USER GUIDES FOR INFORMATION ABOUT
!     PARAMETERS NOT EXPLAINED ABOVE, AND FOR MORE DETAILS ABOUT
!     THE FUNCTION OF THE ROUTINE.
!
!
! **********************************************************************
!
REAL(Kind=WP_Kind)  A, B, CK, CL, COEFF, MU, SMALL, TIMTQL, TLDL, TOPINV
REAL(Kind=WP_Kind)  TOPM, TOTALT, TPRED, TPRED1, TSAVE, TVECOP, W(1), BG, TW, WRR
INTEGER   ADDRSS, CAUSE, CNEGF, CONV, COPT, CPOS, DAFILE, ERRNO
INTEGER   I, ITERNO, IW(1), KFILE, LP, MAXIW, MAXL, MAXW, MXREST
INTEGER   N, NREAD, NUMEIG, NWRITE, PFCONV, PMAX, POPT, RFIRST, RNEW, TCONV
INTEGER   VER, X, TM1, TCSAV, IP1, J, K, TIW, ISTART
INTEGER   WRI, TCONV1, NLEFT, ERRNO1
LOGICAL   F, T, UPDATE, REACHB, USSMXR, USEDB, WRL
COMMON   /STLMAD/ ADDRSS, DAFILE, KFILE, LP, MAXL, NREAD, NWRITE
COMMON   /STLMIN/ A, B, NUMEIG, MAXW, MAXIW
COMMON   /STLMOP/ TLDL, TOPINV, TOPM, TIMTQL, TVECOP, TPRED1, TSAVE, COEFF(4), CK, CL, CONV, PFCONV, UPDATE
COMMON   /STLMPL/ PMAX, POPT, COPT, MXREST
COMMON   /STLMTF/ T, F
COMMON   /STLMTS/ REACHB, USSMXR, USEDB
COMMON   /STLMUI/ MU(3), TOTALT, TPRED, CNEGF, CPOS, ERRNO1,  ITERNO, N, RFIRST, RNEW, TCONV1, X
COMMON   /STLMWR/ WRR(5), WRI(5), WRL(5)
COMMON   /STLMVR/ SMALL, VER
!
!
!     ***********
!     INITIALIZE.
!     ***********
CALL UI
NLEFT = RFIRST
ERRNO = ERRNO1
TCONV = TCONV1
!
!     *********************************************************
!     TERMIN = T, IF WE HAVE REACHED B, THE LAST EIGENVALUE, OR
!                 IF AN ERROR HAS OCCURRED.
!     **********************************************************
TERMIN = (REACHB .AND. USEDB) .OR. CAUSE > 1
IF(REACHB) USEDB = T
IF(.NOT. TERMIN) GOTO 9999
!     ********************************************************
!     BG = LAST SHIFT.
!     IF(REACHED THE LAST EIGENVALUE) BG = B
!     IF(ERROR) BG = THE PREVIOUS SHIFT
!     IF(FIRST ITERATION AND ERROR) BG = A
!     IF(NOT FIRST ITERATION AND CAN NOT COMPUTE A NEW SHIFT)
!           BG = LAST SHIFT
!     ********************************************************
BG = MU(3)
IF(CAUSE == 2) BG = B
IF(CAUSE == 3) BG = MU(2)
IF(ITERNO == 1  .AND. CAUSE == 3) BG = A
IF(ITERNO > 1  .AND.  ERRNO == 203) BG = MU(3)
IF(BG > B) BG = B
IF(BG < A) BG = A
!
!
IF(TCONV <= 0) GOTO 9999
!
!     *********************************************************
!     CONSTRUCT POINTER VECTOR IW, AND MOVE EIGENVALUES TO THE
!     BEGINNING OF W. RESULTS IN W(IW(1)), ..., W(IW(TCONV)) IN
!     ASCENDING ORDER.
!     *********************************************************
DO 5 I = 1, TCONV
  IW(I) = I
  W(I)  = W(LP)
  LP    = LP + 1
5       CONTINUE
!
TM1 = TCONV - 1
IF(TM1 <= 0) GOTO 30
DO 20 I = 1, TM1
  K   = I
  IP1 = I + 1
!
  DO 10 J = IP1, TCONV
    IF(W(J) < W(K)) K = J
10        CONTINUE
!
  IF(.NOT. I /= K) GOTO 20
    TIW   = IW(I)
    IW(I) = IW(K)
    IW(K) = TIW
!
    TW    = W(I)
    W(I)  = W(K)
    W(K)  = TW
20        CONTINUE
!
30    IF(.NOT. ERRNO /= 1701) GOTO 55
!       ***********************************************************
!       IF NOT ERROR DUE TO MAXL TOO SHORT, FIND LARGEST EIGENVALUE
!       SMALLER THAN OR EQUAL TO BG.
!       ***********************************************************
  I = TCONV + 1
  DO 40 J = 1, TCONV
    I = I - 1
    IF(W(I) <= BG) GOTO 50
40        CONTINUE
  I = 0
!
!     **************************************
!     THROW ALL EIGENVALUES GREATER THAN BG.
!     **************************************
50    TCONV = I
!
!
55    IF(.NOT. (MU(1) < A  .AND. TCONV > 0) ) GOTO 9999
  TCSAV = TCONV
!       ****************************************************
!       FIND SMALLEST EIGENVALUE GREATER THAN OR EQUAL TO A.
!       ****************************************************
  DO 60 I = 1, TCONV
    IF(W(I) >= A) GOTO 70
60        CONTINUE
!
  I = TCONV + 1
!
!       *********************************
!       THROW THOSE THAT ARE LESS THAN A.
!       *********************************
70      TCONV = TCONV - (I - 1)
!
!
  ISTART = I
!       ************************************************
!       ADJUST NLEFT, AND MOVE EIGENVALUES AND POINTERS.
!       ************************************************
  NLEFT  = RFIRST + ISTART - 1
  IF(.NOT. (TCONV > 0  .AND.  TCSAV /= TCONV) ) GOTO 9999
    J = 1
    DO 80 I = ISTART, TCSAV
      W(J)  = W(I)
      IW(J) = IW(I)
      J     = J + 1
80          CONTINUE
!
9999  WRI(1) = TCONV
WRI(2) = NLEFT
WRR(1) = BG
!
!
RETURN
END
