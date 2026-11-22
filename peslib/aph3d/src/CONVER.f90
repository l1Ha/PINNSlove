SUBROUTINE CONVER(BETAPI, LAMBDA, POINTR, S, PMAX)
USE Numeric_Kinds_Module
SAVE
!
!   $RCSfile: conver.f,v $ $Revision: 1.3 $
!   $Date: 89/07/28 09:57:06 $
!   $State: Stable $
!
!
! **********************************************************************
!
!     PURPOSE -
!
!         COMPUTE THE NUMBER OF EIGENVALUES IN THE INTERVAL
!         (OLDMU,MU) AND (OLDMU = -INFINITY FOR THE FIRST SHIFT) AND
!         THE NUMBER OF CONVERGED EIGENVALUES IN (MU,INFINITY).
!         COMPUTE CNEG AND CPOS, INITIALIZE POINTR.
!
!     PLEASE SEE THE PROGRAMMERS GUIDE FOR INFORMATION ABOUT
!     PARAMETERS NOT EXPLAINED ABOVE, AND FOR MORE DETAILS ABOUT
!     THE FUNCTION OF THE ROUTINE.
!
!
! **********************************************************************
!
INTEGER   PMAX
REAL(Kind=WP_Kind) ALTMU, BETAPI(1), CK, CL, COEFF, FACTOR, LAMBDA(1),  SRELPR, MU, NEXTMU, OLDMU, S(PMAX, 1), TIMTQL, TLDL
REAL(Kind=WP_Kind) TOLBPI, TOLLDL, TOLPDM, TOLS1I, TOLZBT, TOLZNU, TOPINV, TOPM, TPRED, TSAVE, TVECOP
INTEGER   CNEG, CNEGF, CONV, CPOS, I, ITERNO, LENP, N, OLCPOS, P, PFCONV, POINTR(1), REST, RFIRST, RNEW, ROLD, TCONV
LOGICAL   SAFRST, UPDATE, USEMX, ZERBET
COMMON   /STLMCT/ N, ITERNO, TCONV, CNEG, CPOS, OLCPOS, RNEW, ROLD, REST, P, USEMX, ZERBET
COMMON   /STLMFT/ CNEGF, RFIRST, SAFRST
COMMON   /STLMMU/ MU, OLDMU, NEXTMU, ALTMU
COMMON   /STLMOP/ TLDL, TOPINV, TOPM, TIMTQL, TVECOP, TPRED, TSAVE, COEFF(4), CK, CL, CONV, PFCONV, UPDATE
COMMON   /STLMTL/ SRELPR, TOLBPI, TOLS1I, FACTOR, TOLPDM, TOLZBT, TOLZNU, TOLLDL(3)
!
CNEG = 0
CPOS = 0
LENP = 0
CONV = 0
!
IF(.NOT. ITERNO == 1) GOTO 40
!       ***********************************************************
!       THE FIRST SHIFT IS SPECIAL, SINCE WE DO NOT HAVE ANY OLDMU.
!       ***********************************************************
  DO 30 I = 1, P
    IF(.NOT. ABS(BETAPI(I)) <= TOLBPI) GOTO 30
      CONV = CONV + 1
      IF(.NOT. MU < LAMBDA(I)) GOTO 10
!             ************************************
!             CONVERGED TO THE RIGHT OF THE SHIFT.
!             ************************************
        CPOS = CPOS + 1
        GOTO 20
!
10          IF(.NOT. SAFRST) GOTO 30
!             ***********************************
!             CONVERGED TO THE LEFT OF THE SHIFT.
!             ***********************************
        CNEG = CNEG + 1
!
20          LENP = LENP + 1
!           ******************
!           INITIALIZE POINTR.
!           ******************
      POINTR(LENP) = I
!
30        CONTINUE
!       ***************************
!       LENGTH OF POINTR IS STORED.
!       ***************************
  POINTR(PMAX + 1) = LENP
  GOTO 9999
!
!     ********************************
!     THIS IS NOT THE FIRST ITERATION.
!     ********************************
40    DO 60 I = 1, P
  IF(.NOT. MU < LAMBDA(I)) GOTO 50
    IF(.NOT. ABS(BETAPI(I)) <= TOLBPI) GOTO 60
!           *********
!           AS ABOVE.
!           *********
      CONV = CONV + 1
      CPOS = CPOS + 1
      LENP = LENP + 1
      POINTR(LENP) = I
      GOTO 60
!
50      IF(.NOT. (LAMBDA(I) > OLDMU .AND. REST > 0)) GOTO 60
!         **************************************************************
!         LAMBDA IS IN (OLDMU, MU).
!         IF REST = 0, NOTHING SHOULD CONVERGE, SO WE DO
!         NOT SAVE ANYTHING.
!         IF LAMBDA HAS NOT CONVERGED, POINTR IS GIVEN A NEGATIVE VALUE,
!         BUT CNEG IS NOT UPDATED IN THIS CASE.
!         **************************************************************
    LENP = LENP + 1
    POINTR(LENP) =  - I
    IF(.NOT. ABS(BETAPI(I)) <= TOLBPI) GOTO 60
      POINTR(LENP) = I
      CONV = CONV + 1
      CNEG = CNEG + 1
60      CONTINUE
!
POINTR(PMAX + 1) = LENP
!     *******************************************************
!     IF CNEG =0, UPDATE POINTR, I.E. DELETE NEGATIVE VALUES.
!     *******************************************************
IF(CNEG == 0) CALL CMPRSS(POINTR, PMAX)
!
9999  RETURN
END
