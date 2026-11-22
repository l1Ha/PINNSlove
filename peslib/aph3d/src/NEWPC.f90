SUBROUTINE NEWPC(W, IW)
USE Numeric_Kinds_Module
!
!   $RCSfile: newpc.f,v $ $Revision: 1.3 $
!   $Date: 89/07/28 09:57:26 $
!   $State: Stable $
!
!
! **********************************************************************
!
!     PURPOSE -
!
!         UPDATES POPT, COPT, CK, AND CL. GIVES A PREDICTION ON
!         HOW MUCH TIME THAT WILL BE USED FOR THE NEXT SHIFT.
!
!     PLEASE SEE THE PROGRAMMERS GUIDE FOR INFORMATION ABOUT
!     PARAMETERS NOT EXPLAINED ABOVE, AND FOR MORE DETAILS ABOUT
!     THE FUNCTION OF THE ROUTINE.
!
!
! **********************************************************************
!
REAL(Kind=WP_Kind) TConv, Rold, Wri
REAL(Kind=WP_Kind) CK, CL, COEFF, FCOPT, FPOPT, K, L, SECOND, TEMP
REAL(Kind=WP_Kind) TIMTQL, TLDL, TOPINV, TOPM, TPRED, TSAVE, TVECOP
REAL(Kind=WP_Kind) W(1), WEIGHT, WRR
INTEGER CNEG, CONV, COPT, CPOS, DUMMY, ITERNO, IW(1), LIM, MXREST
INTEGER N, NCOPT, NPOPT, OLCPOS, P, PFCONV, PMAX, POPT,  REST, RNEW, RO
LOGICAL   MEQI, UPDATE, USEMX, WRL, ZERBET, REACHB, USSMXR, USEDB
COMMON   /STLMCT/ N, ITERNO, TCONV, CNEG, CPOS, OLCPOS, RNEW, ROLD,  &
          REST, P, USEMX, ZERBET
COMMON   /STLMOP/ TLDL, TOPINV, TOPM, TIMTQL, TVECOP, TPRED, TSAVE,  &
          COEFF(4), CK, CL, CONV, PFCONV, UPDATE
COMMON   /STLMMI/ MEQI
COMMON   /STLMPL/ PMAX, POPT, COPT, MXREST
COMMON   /STLMTS/ REACHB, USSMXR, USEDB
COMMON   /STLMWR/ WRR(5), WRI(5), WRL(5)
!
IF(.NOT. UPDATE) GOTO 9999
!
WEIGHT = SECOND() - TSAVE
FPOPT = P
FCOPT = CONV
!     ***********************************************************
!     LET T(P, C) = TIME USED (ACCORDING TO OUR FORMULA) FOR ONE
!     SHIFT WHEN LANCZOS IS RUN P STEPS, AND C EIGENPAIRS HAVE
!     CONVERGED.
!     THEN WEIGHT = TIME USED SO FAR FOR THIS SHIFT / T(P, CONV).
!     WEIGHT SHOULD LIE AROUND 1.0, WHICH IS USUALLY DOES.
!     ***********************************************************
WEIGHT = WEIGHT / (COEFF(1) + FPOPT * (COEFF(2) + FCOPT * COEFF(4)+ FPOPT * (TVECOP + FPOPT * TIMTQL)) + FCOPT * COEFF(3))
!
IF(.NOT. CONV > 0) GOTO 10
!       *********************************************
!       COMPUTE NEW K AND L.
!       (P = K * NUMBER OF CONVERGED EIGENPAIRS + L).
!       *********************************************
  L = PFCONV - 1
  TEMP = P
  K = TEMP - L
  TEMP = CONV
  K = K / TEMP
!
!       **************************
!       COMPUTE NEW POPT AND COPT.
!       **************************
  CALL HALF(K, L, NPOPT, NCOPT, LIM)
!
  WRR(4) = WEIGHT
  CALL WRINFO(25, 1, W, IW)
  WRR(1) = L
  WRR(2) = K
!
!       *************************************
!       CHECK THAT THE VALUES ARE REASONABLE.
!       *************************************
  IF(.NOT. (2 * NPOPT >= POPT .AND. NPOPT <= 2 * POPT .AND. &
    2 * NCOPT >= COPT .AND. NCOPT <= 2 * COPT  .AND. LIM >= 4)) GOTO 10
!
!         **************************************
!         YES, THEY ARE. REPLACE THE OLD VALUES.
!         **************************************
    POPT = NPOPT
    COPT = NCOPT
    IF(.NOT. USSMXR) MXREST = COPT
    CK = K
    CL = L
!
10    FPOPT = POPT
FCOPT = COPT
!     ***********************************************************
!     PREDICTED TIME FOR THE NEXT SHIFT = WEIGHT * T(POPT, COPT).
!     ***********************************************************
TPRED = WEIGHT * (COEFF(1) + FPOPT * (COEFF(2) + FCOPT * COEFF(4) + FPOPT * (TVECOP + FPOPT * TIMTQL)) + FCOPT * COEFF(3))
!
CALL WRINFO(25, 2, W, IW)
!
9999  RETURN
END
