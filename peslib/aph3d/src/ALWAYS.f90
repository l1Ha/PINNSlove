LOGICAL FUNCTION ALWAYS(W, IW)
USE Numeric_Kinds_Module
!
!   $RCSfile: always.f,v $ $Revision: 1.3 $
!   $Date: 89/07/28 09:57:03 $
!   $State: Stable $
!
!
! **********************************************************************
!
!     PURPOSE -
!
!         HANDLES THE LANCZOS RUN FOR ONE SHIFT. COMPUTES EIGENVALUES
!         ON BOTH SIDES OF THE SHIFT.
!
!     PLEASE SEE THE PROGRAMMERS GUIDE FOR INFORMATION ABOUT
!     PARAMETERS NOT EXPLAINED ABOVE, AND FOR MORE DETAILS ABOUT
!     THE FUNCTION OF THE ROUTINE.
!
!
! **********************************************************************
!
REAL(Kind=WP_Kind) CK, CL, COEFF, RDUMP, SECOND, ST, TIME, TIMTQL, TLDL, TNORM, TOPINV, TOPM, TPRED, TSAVE, TVECOP, W(1), WRR
INTEGER   ALPHA, BETA, BETAPI, CNEG, CONV, COPT, COUNT, CPOS,  DUMMY, ERRNO, IDUMP, ITERNO, IW(1), LAMBDA, MV, MXNEG,           &
                MXNEW, MXOLD, MXREST, MXRST, N, NBADMU, NIL, NMXRST, NOR, NU, OLCPOS, P, PFCONV, PMAX, POINTR, POPT, REST,      &
                RNEW, ROLD, S, SCPX, SCR, SOLCPX, TCONV, V, WRI, X
LOGICAL   CONTIN, CPALSO, F, FINAL, INITLA, LANCZO, MEQI, POSDOT, POSNU, PREPSV, RANDVC, SAVEXL, T, TRIDIG, UPDATE, USEMX,  &
                VISITL, WRL, ZERBET, REACHB, USSMXR, USEDB
COMMON   /STLMCT/ N, ITERNO, TCONV, CNEG, CPOS, OLCPOS, RNEW, ROLD, REST, P, USEMX, ZERBET
COMMON   /STLMER/ RDUMP, ERRNO, IDUMP(2)
COMMON   /STLMID/ NIL, MV, V, MXNEW, MXOLD, MXRST, SCPX, SOLCPX, X
COMMON   /STLMMI/ MEQI
COMMON   /STLMOP/ TLDL, TOPINV, TOPM, TIMTQL, TVECOP, TPRED, TSAVE, COEFF(4), CK, CL, CONV, PFCONV, UPDATE
COMMON   /STLMPL/ PMAX, POPT, COPT, MXREST
COMMON   /STLMPV/ ALPHA, BETA, BETAPI, LAMBDA, NU, POINTR, S, SCR
COMMON   /STLMST/ TIME(24), COUNT(24), NBADMU, NMXRST, DUMMY
COMMON   /STLMTF/ T, F
COMMON   /STLMTS/ REACHB, USSMXR, USEDB
COMMON   /STLMWR/ WRR(5), WRI(5), WRL(5)
!
DATA      NOR  / 2 /
!
ST = SECOND()
COUNT(NOR) = COUNT(NOR) + 1
ALWAYS = T
!
VISITL = F
!     ******************************************************
!     GENERATE A RANDOMVECTOR AND ORTHOGONALIZE IT AGAINST
!     THE OLCPOS-VECTORS, PRODUCING A STARTINGVECTOR (LET US
!     CALL IT SV).
!     ******************************************************
10    IF(.NOT. RANDVC(V + 1, W, IW, .NOT. USEMX .AND. .NOT. MEQI))GOTO 8888
!
IF(.NOT. PREPSV(SOLCPX, MXOLD, OLCPOS, W, IW)) GOTO 8888
IF(.NOT. INITLA(CONTIN, FINAL, POSDOT, POSNU, W, IW)) GOTO 8888
!     ***********************
!     IS SV(T)*M*SV POSITIVE.
!     ***********************
IF(POSDOT) GOTO 30
  IF(.NOT. VISITL) GOTO 20
    CALL ERROR(NOR, 1)
    GOTO 8888
!
20      VISITL = T
  GOTO 10
!
!     *******************************
!     RUN LANCZOS. CHECK CONVERGENCE.
!     *******************************
30    IF(.NOT. (CONTIN .AND. (.NOT. POSNU .OR.(CNEG == 0 .AND. REST > 0)))) GOTO 40
!
  IF(USEDB .AND. CNEG >= REST) GOTO 50
!
  IF(.NOT. LANCZO(W(ALPHA), W(BETA), TNORM, W, IW)) GOTO 8888
  IF(.NOT. TRIDIG(CONTIN, FINAL, POSNU, T, W, IW)) GOTO 8888
  GOTO 30
!
40    IF(.NOT. (CONTIN .AND. P < POPT .AND.(CNEG < REST .OR. CNEG + CPOS < COPT))) GOTO 45
!
  IF(USEDB .AND. CNEG >= REST) GOTO 50
!
  IF(.NOT. LANCZO(W(ALPHA), W(BETA), TNORM, W, IW)) GOTO 8888
  IF(.NOT. TRIDIG(CONTIN, FINAL, POSNU, T, W, IW)) GOTO 8888
  GOTO 40
!
45    IF(.NOT. (CONTIN .AND. CNEG + CPOS == 0)) GOTO 50
!
  IF(USEDB .AND. CNEG >= REST) GOTO 50
!
  IF(.NOT. LANCZO(W(ALPHA), W(BETA), TNORM, W, IW)) GOTO 8888
  IF(.NOT. TRIDIG(CONTIN, FINAL, POSNU, T, W, IW)) GOTO 8888
  GOTO 45
!
!     **********************
!     NO MORE LANCZOS STEPS.
!     **********************
50    FINAL = T
IF(.NOT. TRIDIG(CONTIN, FINAL, POSNU, T, W, IW)) GOTO 8888
CALL FREEID(MV + P + 1)
!
IF(ITERNO > 1) REST = REST - CNEG
WRR(1) = TNORM
CALL WRINFO(NOR, 1, W, IW)
IF(.NOT. REST < 0) GOTO 60
!       ****************************
!       SHOULD NEVER HAPPEN, BUT ...
!       ****************************
  CALL ERROR(NOR, 2)
  GOTO 8888
!
60    IF(.NOT. REST == 0) GOTO 80
!       ***************************
!       FINISHED WITH THE INTERVAL.
!       ***************************
  IF(.NOT. POSNU) GOTO 70
    CPALSO = T
    MXNEG = NIL
    SCPX = X + TCONV + CNEG
!         *****************************
!         COMPUTE AND STORE EIGENPAIRS.
!         *****************************
    IF(.NOT. SAVEXL(W(BETAPI), W(LAMBDA), W(NU), IW(POINTR), W(S), MXNEG, PMAX, CPALSO, W, IW)) GOTO 8888
!         ******************
!         NEW POPT AND COPT.
!         ******************
    CALL NEWPC(W, IW)
    GOTO 9999
!
!       *****************************************
!       WE DO NOT HAVE A NEW SHIFT, BUT REST = 0.
!       COMPUTE AND STORE EIGENPAIRS.
!       *****************************************
70      CPALSO = F
  MXNEG = NIL
  IF(.NOT. SAVEXL(W(BETAPI), W(LAMBDA), W(NU), IW(POINTR), W(S), MXNEG, PMAX, CPALSO, W, IW)) GOTO 8888
  CALL ERROR(NOR, 3)
  GOTO 8888
!
!     *************************************************
!     WE DO NOT HAVE A NEW SHIFT, AND REST IS POSITIVE.
!     *************************************************
80    IF(POSNU) GOTO 9999
  CALL ERROR(NOR, 4)
!
8888  ALWAYS = F
CALL ERROR(NOR, NOR)
!
9999  TIME(NOR) = TIME(NOR) + SECOND() - ST
RETURN
END FUNCTION ALWAYS
