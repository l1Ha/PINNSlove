SUBROUTINE HALF(K, L, TP, TC, LIM)
USE Numeric_Kinds_Module
USE FileUnits_Module
SAVE
!   $RCSfile: half.f,v $ $Revision: 1.3 $
!   $Date: 89/07/28 09:57:12 $
!   $State: Stable $
!
!
! **********************************************************************
!
!     PURPOSE -
!
!         COMPUTES COEFFICIENTS IN THE POLYNOMIAL WHOSE ZERO GIVES
!         COPT. COMPUTES THE ZERO AND POPT.
!
!
!     INPUT PARAMETERS -
!
!         K,L =  WE ASSUME THAT P=K*C+L, WHERE C = NUMBER OF CONVERGED
!                EIGENPAIRS FOR P LANCZOS STEPS.
!
!
!     OUTPUT PARAMETERS -
!
!         TP  = TEMPORARY POPT.
!         TC  = TEMPORARY COPT.
!         LIM = 6 NORMAL TERMINATION.
!               1 K OR L NEGATIVE.
!               2 CAN NOT FIND COPT.
!               3 TP OR TL IS ABSURD.
!               4 PMAX IS TOO SMALL COMPARED WITH TP.
!               5 TP IS GREATER THAN N.
!
!     PLEASE SEE THE PROGRAMMERS GUIDE FOR INFORMATION ABOUT
!     PARAMETERS NOT EXPLAINED ABOVE, AND FOR MORE DETAILS ABOUT
!     THE FUNCTION OF THE ROUTINE.
!
!
! **********************************************************************
!
REAL(Kind=WP_Kind) A, B, C, CK, CL, COEFF, FM, FN, FP, K, L, MEAN, NEG, POS
REAL(Kind=WP_Kind) TEMP, TIMTQL, TLDL, TOPINV, TOPM, TPRED, TSAVE, TVECOP, WRR
INTEGER   CNEG, CONV, COPT, CPOS, ITERNO, LIM, LOOP, MXREST, N,  OLCPOS, P, PFCONV, PMAX, POPT, REST, RNEW, ROLD, TC, TCONV, TP, WRI
LOGICAL   MEQI, UPDATE, USEMX, WRL, ZERBET
COMMON   /STLMCT/ N, ITERNO, TCONV, CNEG, CPOS, OLCPOS, RNEW, ROLD, REST, P, USEMX, ZERBET
COMMON   /STLMOP/ TLDL, TOPINV, TOPM, TIMTQL, TVECOP, TPRED, TSAVE, COEFF(4), CK, CL, CONV, PFCONV, UPDATE
COMMON   /STLMMI/ MEQI
COMMON   /STLMPL/ PMAX, POPT, COPT, MXREST
COMMON   /STLMWR/ WRR(5), WRI(5), WRL(5)
!
LIM = 6
IF(.NOT. (K <= 0.0D0 .OR. L <= 0.0D0)) GOTO 10
!       **************************
!       SHOULD NOT HAPPEN, BUT ...
!       **************************
  LIM = 1
  GOTO 9999
!
!     **************************************************
!     COMPUTE COEFFICIENTS IN TIME OPTIMIZATION FORMULA.
!     **************************************************
10    A = 2.0D0 * TIMTQL * K**3
B = K * COEFF(4) + TVECOP * K * K + 3.0D0 * TIMTQL * K * K * L
C = COEFF(1) + L * (COEFF(2) + L * (TVECOP + TIMTQL * L))
WRR(1) = A
WRR(2) = B
WRR(3) = C
WRITE(Out_Unit,*)'c=',c,' b=',b
IF(B<=0.d0)b=1.d0  !tempmod
POS = SQRT(C / B)
FP = A * POS**3
NEG = 0.0D0
FN =  - C
LOOP = 0
!
!     ********************************************************
!     SEEK AN APPROXIMATION TO A ZERO TO
!     A * X**3 + B * X**2 - C = 0, BY SUCCESSIVE BISECTIONS.
!     LOOP IS USED TO AVOID AN INFINITE LOOP IN CASE SOMETHING
!     SHOULD GO WRONG. (LOOP = 5 IS NORMAL).
!     ********************************************************
20    IF(.NOT. (POS - NEG > 4.5D-1 / K .AND. LOOP <= 100)) GOTO 50
  MEAN = (NEG + POS) * 0.5D0
  FM = MEAN * MEAN * (A * MEAN + B) - C
  IF(FM == 0.0D0) GOTO 70
  IF(.NOT. FM > 0.0D0) GOTO 30
    POS = MEAN
  GOTO 40
30        NEG = MEAN
40      LOOP = LOOP + 1
  GOTO 20
!
50    IF(.NOT. LOOP > 100) GOTO 60
  LIM = 2
  GOTO 9999
!
60    MEAN = (NEG + POS) * 0.5D0
!     ********************************
!     COMPUTE TEMPORARY POPT AND COPT.
!     ********************************
70    TC = MEAN + 0.5D0
TEMP = TC
TP = K * TEMP + L + 0.5D0
!     ********************************************************
!     EXPERIENCE HAS SHOWN THAT A TP < 20 CAN CAUSE TROUBLE
!     WITH THE SHIFTING.
!     ********************************************************
IF(.NOT. TP < 20) GOTO 75
  TP = 20
  TEMP = TP
  TC = MAX(1.5D0, (TEMP - L) / K +0.5D0)
75    WRI(1) = TC
WRI(2) = TP
WRI(3) = LOOP
!
IF(.NOT. (TC <= 0 .OR. TP <= 0)) GOTO 80
!       ****************************
!       SHOULD NEVER HAPPEN, BUT ...
!       ****************************
  LIM = 3
  GOTO 9999
!
80    IF(.NOT. TP > PMAX) GOTO 90
!       ******************************************
!       WE CAN NOT TAKE SO MANY STEPS DUE TO PMAX.
!       ******************************************
  LIM = 4
  TP = PMAX
  TEMP = TP
  TC = MAX(1.5D0, (TEMP - L) / K + 0.5D0)
!
90    IF(.NOT. TP > N) GOTO 9999
!       ***************************************
!       THE PROBLEM IS VERY SMALL, LIMITS POPT.
!       ***************************************
  LIM = 5
  TP = N
  TEMP = TP
  TC = MAX(1.5D0, (TEMP - L) / K + 0.5D0)
!
9999  WRI(4) = LIM
RETURN
END
