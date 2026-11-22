SUBROUTINE DELDUP(NU, POINTR, S, MAXNU, MINNU, NUMDEL, PMAX,  FINAL)
USE Numeric_Kinds_Module
SAVE
!
!   $RCSfile: deldup.f,v $ $Revision: 1.3 $
!   $Date: 89/07/28 09:57:08 $
!   $State: Stable $
!
!
! **********************************************************************
!
!     PURPOSE -
!
!         DELETES DUPLICATE EIGENPAIRS TO THE LEFT OF THE SHIFT. THERE
!         ARE TWO CASES, VIZ.
!
!         T HAS NUMERICALLY MULTIPLE EIGENVALUES. LIKELY TO HAPPEN IF
!         WE ARE CLOSE TO AN EIGENVALUE OF (HIGH) MULTIPLICITY.
!
!         THE EIGENVECTOR HAS A SMALL COMPONENT IN THE DIRECTION OF
!         THE STARTINGVECTOR. LIKELY TO HAPPEN IF MU IS CLOSE TO THE
!         CORRESPONDING EIGENVALUE AND THE PAIR IS ALREADY COMPUTED.
!         WE CHECK ON S(1,I).
!
!
!     INPUT PARAMETERS -
!
!         MAXNU  = MAX(NU(I)), I=1, ..., P.
!         MINNU  = MIN(NU(I)), I=1, ..., P.
!         NUMDEL = NUMBER OF DELETED VECTORS IN CASE 1.
!         FINAL  = TRUE IF THE LAST LANCZOS STEP HAS BEEN TAKEN.
!
!     PLEASE SEE THE PROGRAMMERS GUIDE FOR INFORMATION ABOUT
!     PARAMETERS NOT EXPLAINED ABOVE, AND FOR MORE DETAILS ABOUT
!     THE FUNCTION OF THE ROUTINE.
!
!
! **********************************************************************
!
INTEGER   PMAX
REAL(Kind=WP_Kind)  ADNU, FACTOR, SRELPR, MAXNU, MINNU, NU(1), NU1, NU2, S1, S2, S(PMAX, 1), TOL, TOLBPI, TOLDEL, TOLLDL, TOLPDM,  &
                TOLS1I, TOLZBT, TOLZNU, TSMAX, WRR
INTEGER   CNEG, CPOS, I, ITERNO, LEN, N, NUMDEL, OLCPOS, P,  POINT, POINTR(1), REST, RNEW, ROLD, S1IDEL, TCONV, WRI
LOGICAL   FINAL, USEMX, WRL, ZERBET
COMMON   /STLMCT/ N, ITERNO, TCONV, CNEG, CPOS, OLCPOS, RNEW, ROLD, REST, P, USEMX, ZERBET
COMMON   /STLMTL/ SRELPR, TOLBPI, TOLS1I, FACTOR, TOLPDM, TOLZBT, TOLZNU, TOLLDL(3)
COMMON   /STLMWR/ WRR(5), WRI(5), WRL(5)
!
!     ************************************************
!     NUMDEL = NUMBER OF DELETED PAIRS, FIRST REASON.
!     S1IDEL =                        , SECOND REASON.
!     ************************************************
NUMDEL = 0
S1IDEL = 0
IF(CNEG == 0 .OR. ITERNO == 1) GOTO 9999
!     **********************
!     COMPUTE THE TOLERANCE.
!     **********************
TOLDEL = (MAXNU - MINNU) * FACTOR
TOL = TOLDEL * 1.0D-4
IF(MAXNU == MINNU) TOLDEL = 1.0D0
I = 1
!     ***********************************
!     IABS, SINCE POINTR MAY BE NEGATIVE.
!     ***********************************
POINT = IABS(POINTR(I))
NU2 = NU(POINT)
!     *****************************************************************
!     S1 AND S2 ARE THE CORRESPONDING TOP ELEMENTS OF THE EIGENVECTORS.
!     *****************************************************************
S2 = ABS(S(1, POINT))
!
!     *******************************************
!     LEN = THE NUMBER OF LAMBDAS IN (OLDMU, MU).
!     *******************************************
LEN = POINTR(PMAX + 1) - CPOS
10    IF(.NOT. I < LEN) GOTO 40
  NU1 = NU2
  S1 = S2
  I = I + 1
  POINT = IABS(POINTR(I))
  NU2 = NU(POINT)
  S2 = ABS(S(1, POINT))
!
  ADNU = ABS(NU1 - NU2)
!       ********************************************************
!       THE CONTROL CONSISTS OF TWO STEPS.
!         IS ADNU SMALL ENOUGH (WE MUST HAVE THIS CONTROL SINCE
!         THE SECOND MAY BE TRUE IF ADNU IS BIG).
!
!         CHECK ON ADNU*(THE QUOTIENT BETWEEN THE TOP ELEMENTS).
!
!       TO AVOID DIVISION WITH ZERO, WE DO NOT FORM THE
!       QUOTIENT EXPLICITLY.
!       ********************************************************
  IF(.NOT. ADNU <= TOLDEL) GOTO 10
!
    TSMAX = MAX(S1, S2)
    IF(.NOT. ADNU * MIN(S1, S2) <= TSMAX * TOL) GOTO 10
      IF(TSMAX == 0.0D0) GOTO 10
!
!           ******************************************************
!           CHECK SO THAT WE DO NOT HAVE TWO NOT CONVERGED VALUES.
!           ******************************************************
      IF(POINTR(I - 1) < 0 .AND. POINTR(I) < 0) GOTO 10
!
!           ******************
!           MAKE THE ROTATION.
!           ******************
      CALL ROTATE(S, POINTR(I - 1), FINAL, PMAX, NUMDEL)
!
      I = I + 1
      IF(I >= LEN) GOTO 40
      POINT = IABS(POINTR(I))
      NU2 = NU(POINT)
      S2 = ABS(S(1, POINT))
      GOTO 10
!
!     ********************************
!     DELETE THOSE WITH SMALL S(1, I).
!     ********************************
40    DO 50 I = 1, LEN
  POINT = POINTR(I)
  IF(POINT <= 0) GOTO 50
  IF(.NOT. ABS(S(1, POINT)) <= TOLS1I) GOTO 50
    POINTR(I) = 0
    S1IDEL = S1IDEL + 1
50      CONTINUE
!
!     **************
!     UPDATE POINTR.
!     **************
IF(FINAL) CALL CMPRSS(POINTR, PMAX)
!
!
9999  WRI(2) = S1IDEL
!     ************
!     ADJUST CNEG.
!     ************
CNEG = CNEG - NUMDEL - S1IDEL
!
RETURN
END
