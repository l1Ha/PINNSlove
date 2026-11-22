LOGICAL FUNCTION TRIDIG(CONTIN, FINAL, POSNU, UPDMU, W, IW)
USE Numeric_Kinds_Module
!
!   $RCSfile: tridig.f,v $ $Revision: 1.3 $
!   $Date: 89/07/28 09:57:49 $
!   $State: Stable $
!
!
! **********************************************************************
!
!     PURPOSE -
!
!         THIS ROUTINE HANDLES ALL COMPUTATIONS CONNECTED WITH THE
!         EIGENSYSTEM OF T (THE TRIDIAGONAL P X P-MATRIX).
!
!
!     INPUT PARAMETERS -
!
!         FINAL = TRUE IF WE HAVE TAKEN THE LAST LANCZOS STEP, FALSE
!                 OTHERWISE.
!         UPDMU = TRUE, IF NEXTMU SHOULD BE COMPUTED, FALSE
!                 OTHERWISE.
!
!
!     OUTPUT PARAMETERS -
!
!         CONTIN = TRUE IF WE CAN TAKE AN OTHER STEP, I.E. IF THERE IS
!                  ROOM LEFT (P < PMAX), VECTORS HAVE NOT BEGUN TO
!                  SPLIT (NUMDEL == 0), AND BETA(P) IS LARGE ENOUGH
!                  (.NOT. ZERBET). FALSE OTHERWISE.
!         POSNU  = TRUE IF WE HAVE A POSITIVE NU EIGENVALUE OR IF MU IS
!                  GREATER THAN LAMBDA(N), OR IF WE HAVE USED
!                  B AS SHIFT. FALSE OTHERWISE.
!
!     PLEASE SEE THE PROGRAMMERS GUIDE FOR INFORMATION ABOUT
!     PARAMETERS NOT EXPLAINED ABOVE, AND FOR MORE DETAILS ABOUT
!     THE FUNCTION OF THE ROUTINE.
!
!
! **********************************************************************
!
REAL(Kind=WP_Kind) CK, CL, COEFF, MAXNU, MINNU, SECOND, ST, TIME, TIMTQL,  TLDL, TOPINV, TOPM, TPRED, TSAVE, TVECOP, W(1), WRR
INTEGER   ALPHA, BETA, BETAPI, CNEG, CNEGF, CONV, COPT, COUNT, CPOS, DUMMY, ITERNO, IW(1), LAMBDA, MXREST, N, NBADMU,         &
          NMXRST, NOR, NU, NUMDEL, OLCPOS, P, PFCONV, PMAX, POINTR, POPT, REST, RFIRST, RNEW, ROLD, S, SCR, TCONV, WRI
LOGICAL   CONTIN, F, FINAL, IMTQL2, POSNU, SAFRST, T, UPDATE,  UPDMU, USEMX, WRL, ZERBET, REACHB, USSMXR, USEDB, LGEB
COMMON   /STLMCT/ N, ITERNO, TCONV, CNEG, CPOS, OLCPOS, RNEW, ROLD, REST, P, USEMX, ZERBET
COMMON   /STLMFT/ CNEGF, RFIRST, SAFRST
COMMON   /STLMOP/ TLDL, TOPINV, TOPM, TIMTQL, TVECOP, TPRED, TSAVE, COEFF(4), CK, CL, CONV, PFCONV, UPDATE
COMMON   /STLMPL/ PMAX, POPT, COPT, MXREST
COMMON   /STLMPV/ ALPHA, BETA, BETAPI, LAMBDA, NU, POINTR, S, SCR
COMMON   /STLMST/ TIME(24), COUNT(24), NBADMU, NMXRST, DUMMY
COMMON   /STLMTS/ REACHB, USSMXR, USEDB
COMMON   /STLMTF/ T, F
COMMON   /STLMWR/ WRR(5), WRI(5), WRL(5)
!
DATA      NOR  / 23 /
!
ST = SECOND()
COUNT(NOR) = COUNT(NOR) + 1
TRIDIG = T
!
!     *******************************
!     SOLVE TRIDIAGONAL EIGENPROBLEM.
!     *******************************
IF(.NOT. IMTQL2(W(ALPHA), W(BETA), W(NU), W(SCR), W(S), P, PMAX, FINAL)) GOTO 8888
!
!     **************************
!     COMPUTE LAMBDA AND BETAPI.
!     **************************
CALL COMPL(W(BETA), W(BETAPI), W(NU), W(LAMBDA), W(S), MAXNU, MINNU, PMAX, LGEB)
!
!     ******************
!     CHECK CONVERGENCE.
!     ******************
CALL CONVER(W(BETAPI), W(LAMBDA), IW(POINTR), W(S), PMAX)
!
!     ************
!     SORT POINTR.
!     ************
CALL SORTP(IW(POINTR), W(LAMBDA), PMAX + 1)
!
IF(FINAL) CALL WRINFO(NOR, 1, W, IW)
!
!     ***********************
!     DELETE DUPLICATE PAIRS.
!     ***********************
CALL DELDUP(W(NU), IW(POINTR), W(S), MAXNU, MINNU, NUMDEL, PMAX,  FINAL)
!
CONV = CONV - NUMDEL
IF(CONV > 0 .AND. PFCONV == 0) PFCONV = P
!     ******************************************************
!     NOTE RNEW == N. WE CAN NOT EXPECT TO GET A NEW SHIFT
!     IF WE ARE TO THE RIGHT OF LAMBDA(N).
!     ******************************************************
POSNU = MAXNU > 0.0D0 .OR. RNEW == N .OR. USEDB
!
!     **********************************************
!     CHECK IF WE SHOULD TAKE AN OTHER LANCZOS STEP.
!     **********************************************
CONTIN = P < N .AND. P < PMAX .AND. NUMDEL == 0 .AND. .NOT. ZERBET
IF(.NOT. FINAL) GOTO 9999
!       ***************************************************
!       NO MORE LANCZOS STEPS. RNEW == N, SEE NOTE ABOVE.
!       COMPUTE THE NEW SHIFT.
!       ***************************************************
  IF(RNEW == N) CPOS = 0
  IF(UPDMU .AND. MAXNU > 0.0D0 .AND. RNEW < N   .AND. (.NOT. USEDB) )  CALL NEWMU(W(LAMBDA), IW(POINTR), MAXNU, LGEB)
!
  IF(ITERNO == 1) CNEGF = CNEG
!
  WRI(1) = NUMDEL
  WRL(1) = POSNU
  WRR(1) = MINNU
  WRR(2) = MAXNU
  CALL WRINFO(NOR, 2, W, IW)
!
  GOTO 9999
!
8888  TRIDIG = F
CALL ERROR(NOR, NOR)
!
9999  TIME(NOR) = TIME(NOR) + SECOND() - ST
RETURN
END
