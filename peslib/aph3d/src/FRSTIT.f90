LOGICAL FUNCTION FRSTIT(W, IW)
USE Numeric_Kinds_Module
!
!   $RCSfile: frstit.f,v $ $Revision: 1.3 $
!   $Date: 89/07/28 09:57:11 $
!   $State: Stable $
!
!
! **********************************************************************
!
!     PURPOSE -
!
!         DETERMINES WHAT METHOD SHOULD BE USED FOR THE
!         ORTHOGONALIZATION OF THE STARTINGVECTOR. INITIALIZES POPT,
!         COPT, AND CONSTANTS USED IN LATER UPDATES. GIVES VALUES TO
!         /STLMID/.
!
!     PLEASE SEE THE PROGRAMMERS GUIDE FOR INFORMATION ABOUT
!     PARAMETERS NOT EXPLAINED ABOVE, AND FOR MORE DETAILS ABOUT
!     THE FUNCTION OF THE ROUTINE.
!
!
! **********************************************************************
!
REAL(Kind=WP_Kind) CK, CL, COEFF, DOT, SECOND, ST, TCK, TCL, TEMP, TIME, TIMTQL, TLDL, Point4
REAL(Kind=WP_Kind) TOPINV, TOPM, TPRED, TSAVE, TVECOP, W(1)
INTEGER  ADR, ALPHA, BETA, BETAPI, CNEG, CNEGF, CONV, COPT, countit, CPOS, DUMMY, FREE
INTEGER  I, ITERNO, IW(1), L, LAMBDA, LIM, MV, MXNEW, MXOLD, MXREST, MXRST, N, NBADMU
INTEGER  NIL, NMXRST, NOACTN, NOR, NU, OLCPOS, P, PFCONV, PMAX, POINTR, POPT, REST
INTEGER  RFIRST, RNEW, ROLD, S, saveit, SAVFRE, SCPX, SCR, SOLCPX, TCONV, TCOPT, TPOPT, V, X
LOGICAL  F, FINAL, IMTQL2, MEQI, MULVEC, OPINV, OPM, RANDVC, SAFRST, SCPROD, SUBVEC
LOGICAL  T, UPDATE, USEID, USEMX, USEMX1, REACHB, USSMXR, USEDB
COMMON   /STLMAC/ NOACTN, FREE, saveit, SAVFRE
COMMON   /STLMCT/ N, ITERNO, TCONV, CNEG, CPOS, OLCPOS, RNEW,  ROLD, REST, P, USEMX, USEMX1
COMMON   /STLMFT/ CNEGF, RFIRST, SAFRST
COMMON   /STLMID/ NIL, MV, V, MXNEW, MXOLD, MXRST, SCPX, SOLCPX, X
COMMON   /STLMMI/ MEQI
COMMON   /STLMOP/ TLDL, TOPINV, TOPM, TIMTQL, TVECOP, TPRED, TSAVE, COEFF(4), CK, CL, CONV, PFCONV, UPDATE
COMMON   /STLMPL/ PMAX, POPT, COPT, MXREST
COMMON   /STLMPV/ ALPHA, BETA, BETAPI, LAMBDA, NU, POINTR, S, SCR
COMMON   /STLMST/ TIME(24), countit(24), NBADMU, NMXRST, DUMMY
COMMON   /STLMTS/ REACHB, USSMXR, USEDB
COMMON   /STLMTF/ T, F
!
DATA      NOR  / 5 /
!
ST = SECOND()
countit(NOR) = countit(NOR) + 1
FRSTIT = T
!
RFIRST = RNEW
IF(.NOT. (RFIRST == N .AND. .NOT. SAFRST)) GOTO 10
!       *********************************************************
!       THE FIRST SHIFT IS GREATER THAN THE LAST EIGENVALUE, AND
!       WE DO NOT WISH TO SAVE THOSE THAT MIGHT CONVERGE TO THE
!       LEFT OF THE SHIFT. HE MUST IN THE TOMS VERSION.
!       *********************************************************
  CALL ERROR(NOR, 1)
  GOTO 8888
!
10    USEID = T
!     *************************************
!     CHECK IF THE USER HAS GIVEN /STLMID/.
!     HE MAY NOT IN THE TOMS VERSION.
!     *************************************
IF(.NOT. V == X) GOTO 15
  USEID = F
  MV = 0
  V = PMAX + 1
  X = V + 2
  IF(MEQI) V = MV
!
!     *********************************************
!     TLDL = TIME FOR THE FIRST LDLT-DECOMPOSITION.
!     *********************************************
15    TLDL = TIME(3) - TLDL
!     ****************************************************************
!     IF( DEFAULT MX OR UPDATE OF POPT, COPT )THEN ... ELSE GOTO 100.
!     ****************************************************************
IF(.NOT. ((USEMX .AND. .NOT. USEMX1) .OR. UPDATE)) GOTO 100
  IF(.NOT. RANDVC(V + 1, W, IW, F)) GOTO 8888
  IF(.NOT. TOPM <= 0.0D0) GOTO 20
    TOPM = 0.0D0
    IF(MEQI) GOTO 20
!         ************************************************************
!         TOPM = TIME FOR ONE Y = M * X. THE CONTROL ABOVE IS
!         MADE TO SEE IF THE USER HAS SUPPLIED A VALUE ON TOPM OR NOT.
!         HE MAY NOT SUPPLY VALUES IN THE TOMS VERSION.
!         ************************************************************
    TOPM = SECOND()
    IF(.NOT. OPM(MV + 1, V + 1, W, IW)) GOTO 8888
    TOPM = SECOND() - TOPM
!
20      IF(.NOT. TVECOP <= 0.0D0) GOTO 40
    TEMP = 3 * N
!         *********************************************
!         HOW MANY LOOPS. ADAPTED FOR A MACHINE WHERE A
!         MULTIPLICATION TAKES ABOUT 1E-5 SECONDS.
!         *********************************************
    L = MAX(Three/Two, One/Two / (TEMP * 1.0D-5))
    IF(L > 300) L = 300
    TVECOP = SECOND()
!         ***************************
!         TIME THE VECTOR OPERATIONS.
!         ***************************
    ADR = saveit
    IF(MEQI) ADR = NOACTN
    Point4=Two/Five
    DO 30 I = 1, L
      IF(.NOT. MULVEC(MV + 1, V + 1, Point4, W, IW, ADR, 1)) GOTO 8888
      IF(.NOT. SUBVEC(V + 1, MV + 1, Point4, W, IW, saveit, 1)) GOTO 8888
      IF(.NOT. SCPROD(MV + 1, V + 1, DOT, W, IW, saveit, 1)) GOTO 8888
30          CONTINUE
    TEMP = 3 * L
!         ****************************************
!         COMPUTE THE MEAN VALUE FOR N OPERATIONS.
!         ****************************************
    TVECOP = (SECOND() - TVECOP) / TEMP
!
!       ****************************************************
!       IF ONLY MX-DECISION, WE DO NOT NEED THE TIMES BELOW.
!       ****************************************************
40      IF(.NOT. UPDATE) GOTO 90
    IF(.NOT. TOPINV <= 0.0D0) GOTO 50
!           ********************************************
!           TOPINV = TIME FOR ONE X =(LDL(T)**(-1)) * B.
!           ********************************************
      TOPINV = SECOND()
      IF(.NOT. OPINV(MV + 1, V + 1, W, IW)) GOTO 8888
      TOPINV = SECOND() - TOPINV
!
50        IF(.NOT. TIMTQL <= 0.0D0) GOTO 70
!           *********************************************
!           SIMULATE THE TIME FOR ONE COMPUTATION OF THE
!           EIGENPAIRS OF THE TRIDIAGONAL MATRIX T.
!           *********************************************
      L = MIN0(PMAX, 20)
      FINAL = F
      TIMTQL = SECOND()
      DO 60 I = 1, L
        ADR = ALPHA + I - 1
        W(ADR) = I
        ADR = BETA + I - 1
        W(ADR) = 1.0D0
        IF(I == L) FINAL = T
        IF(.NOT. IMTQL2(W(ALPHA), W(BETA), W(NU), W(SCR), W(S), I, PMAX, FINAL)) GOTO 8888
60            CONTINUE
      TIMTQL = SECOND() - TIMTQL
!
    TEMP = L**3
!         ********************************************
!         COMPUTE C WHERE TIMTQL = C * (DIM. OF T)**3.
!         COMPUTE COEFFICIENTS USED IN NEWPC AND HALF.
!         ********************************************
    TIMTQL = TIMTQL / TEMP
70        COEFF(1) = TLDL + 2.0D0 * (TOPM + TVECOP)
    IF(MEQI) COEFF(1) = COEFF(1) - TVECOP
    COEFF(2) = TOPINV + 7.0D0 * TVECOP + TOPM
    IF(MEQI) COEFF(2) = COEFF(2) - TVECOP
    COEFF(3) = 2.0D0 * TVECOP
    COEFF(4) = TVECOP
!
!         ******************************************
!         COMPUTE POPT AND COPT, AND CHECK IF OK.
!         THE TEMPORARY VARIABLES ARE TO PASS PFORT.
!         ******************************************
    TCK = CK
    TCL = CL
    TCOPT = COPT
    TPOPT = POPT
    CALL HALF(TCK, TCL, TPOPT, TCOPT, LIM)
    COPT = TCOPT
    POPT = TPOPT
    IF(.NOT. LIM <= 3) GOTO 80
      POPT = MIN0(N, PMAX)
      TEMP = PMAX - 10
      COPT = TEMP * 0.5D0 + 0.5D0
      IF(COPT < 1) COPT = 1
!
80        IF(.NOT. USSMXR) MXREST = COPT
!
90      CALL FREEID(V + 1)
  CALL FREEID(MV + 1)
!
100   TEMP = POPT + 1
!     ******************************
!     SHOULD WE SAVE THE MV-VECTORS.
!     ******************************
IF(USEMX .AND. .NOT. USEMX1) USEMX = TEMP * TVECOP < TOPM
IF(USEMX .AND. MEQI) USEMX = F
IF(USEID) GOTO 110
!       *************************************************
!       INITIALIZE /STLMID/ IF THE USER HAS NOT DONE IT.
!       HE MAY NOT DO IT IN THE TOMS VERSION.
!       *************************************************
  X = V + PMAX + 1
  IF(.NOT. USEMX) GOTO 110
    MXNEW = V + PMAX + 1
    MXOLD = MXNEW + PMAX
    MXRST = MXOLD + PMAX
    X = MXRST + PMAX
!
110   IF(.NOT. UPDATE) GOTO 9999
!       ******************
!       CORRECT FOR USEMX.
!       ******************
  I = 0
  IF(USEMX) I = 1
  TEMP = I
  COEFF(3) = COEFF(3) + (1.0D0 - TEMP) * TOPM * 0.5D0 + TEMP * TVECOP * 0.5D0
  COEFF(4) = COEFF(4) + TEMP * TVECOP * 0.5D0
  GOTO 9999
!
8888  FRSTIT = F
CALL ERROR(NOR, NOR)
!
9999  IF(FRSTIT) CALL WRINFO(NOR, 1, W, IW)
TIME(NOR) = TIME(NOR) + SECOND() - ST
TSAVE = TSAVE + TIME(NOR)
RETURN
END
