LOGICAL FUNCTION SELDOM(W, IW)
USE Numeric_Kinds_Module
!
!   $RCSfile: seldom.f,v $ $Revision: 1.3 $
!   $Date: 89/07/28 09:57:40 $
!   $State: Stable $
!
!
! **********************************************************************
!
!     PURPOSE -
!
!         HANDLES THE CASE WHEN REST IS POSITIVE AFTER ALWAYS, I.E.
!         THERE ARE MISSING EIGENVALUES IN THE INTERVAL. THE ALGORITHM
!         IS AS FOLLOWS,
!
!               10 ORTHOGONALIZE THE STARTINGVECTOR AGAINST THE VECTORS
!                  THAT HAVE CONVERGED, AND WHOSE EIGENVALUES LIE
!                  IN (OLDMU,MU) (AND AGAINST THOSE TO THE RIGHT OF
!                  MU THE FIRST TIME).
!
!                  RUN LANCZOS AS USUAL.
!
!                  SAVE THE NEW VECTORS WITH LAMBDA IN THE INTERVAL.
!
!                  IF STILL MISSING, GOTO 10.
!
!     PLEASE SEE THE PROGRAMMERS GUIDE FOR INFORMATION ABOUT
!     PARAMETERS NOT EXPLAINED ABOVE, AND FOR MORE DETAILS ABOUT
!     THE FUNCTION OF THE ROUTINE.
!
!
!
! **********************************************************************
!
REAL(Kind=WP_Kind)  RDUMP, SECOND, ST, TIME, TNORM, W(1), WRR
INTEGER   ALPHA, BETA, BETAPI, CNEG, COPT, COUNT, CPOS, DUMMY, ERRNO, IDUMP, ITERNO, IW(1), LAMBDA, MV, MXNEG, MXNEW,          &
                MXOLD, MXREST, MXRST, N, NBADMU, NIL, NMXRST, NOR, NU, NUMR, OLCPOS, P, PMAX, POINTR, POPT, REST, RNEW, ROLD,  &
                S, SACNEG, SACPOS, SCNEGX, SCPOSX, SCPX, SCR, SOLCPX, SRESTX, TCONV, V, WRI, X
LOGICAL   CONTIN, CPALSO, F, FINAL, INITLA, LANCZO, MEQI, POSDOT, POSNU, PREPSV, RANDVC, SAVEXL, T, TRIDIG, USEMX, VISITL,     &
                WRL, ZERBET
COMMON   /STLMCT/ N, ITERNO, TCONV, CNEG, CPOS, OLCPOS, RNEW,  ROLD, REST, P, USEMX, ZERBET
COMMON   /STLMER/ RDUMP, ERRNO, IDUMP(2)
COMMON   /STLMID/ NIL, MV, V, MXNEW, MXOLD, MXRST, SCPX, SOLCPX, X
COMMON   /STLMMI/ MEQI
COMMON   /STLMPL/ PMAX, POPT, COPT, MXREST
COMMON   /STLMPV/ ALPHA, BETA, BETAPI, LAMBDA, NU, POINTR, S, SCR
COMMON   /STLMST/ TIME(24), COUNT(24), NBADMU, NMXRST, DUMMY
COMMON   /STLMTF/ T, F
COMMON   /STLMWR/ WRR(5), WRI(5), WRL(5)
!
DATA      NOR  / 19 /
!
ST = SECOND()
COUNT(NOR) = COUNT(NOR) + 1
SELDOM = T
!
!     ************************************************************
!     SET UP POINTER INFORMATION.
!       SCNEGX = START OF CNEG X-VECTORS.
!       SCPOSX =          CPOS
!       SRESTX =          REST
!
!     SAVE CNEG AND CPOS, SINCE THEY WOULD BE DESTROYED OTHERWISE.
!     NUMR = NUMBER OF REST VECTORS.
!     ************************************************************
SCNEGX = X + TCONV
SCPOSX = SCNEGX + CNEG
SRESTX = SCPOSX + CPOS
SACNEG = CNEG
SACPOS = CPOS
NUMR = 0
!
CPALSO = T
!     ********************************************************
!     MXNEG = START OF CNEG MX-VECTORS.
!     NOTE, THERE IS ROOM FOR THE CPOS MX-VECTORS TO THE LEFT.
!     ********************************************************
MXNEG = MXNEW + CPOS
SCPX = SCPOSX
!     *********************************
!     COMPUTE AND STORE THE EIGENPAIRS.
!     *********************************
IF(.NOT. SAVEXL(W(BETAPI), W(LAMBDA), W(NU), IW(POINTR), W(S), MXNEG, PMAX, CPALSO, W, IW)) GOTO 8888
!
VISITL = F
!
!     *************************************************
!     COMPUTE RANDOMVECTOR AND ORTHOGONALIZE IT AGAINST
!       CPOS
!       OLD CPOS
!       CNEG
!       REST X-VECTORS.
!     *************************************************
10    IF(.NOT. RANDVC(V + 1, W, IW, .NOT. USEMX .AND. .NOT. MEQI)) GOTO 8888
!
IF(.NOT. PREPSV(SCPOSX, MXNEW, SACPOS, W, IW)) GOTO 8888
IF(.NOT. PREPSV(SOLCPX, MXOLD, OLCPOS, W, IW)) GOTO 8888
IF(.NOT. PREPSV(SCNEGX, MXNEW + SACPOS, SACNEG, W, IW)) GOTO 8888
IF(.NOT. PREPSV(SRESTX, MXRST, NUMR, W, IW)) GOTO 8888
IF(.NOT. INITLA(CONTIN, FINAL, POSDOT, POSNU, W, IW)) GOTO 8888
!
IF(POSDOT) GOTO 30
!       ************************************
!       M INDEFINITE, OR BAD STARTINGVECTOR.
!       ************************************
  IF(.NOT. VISITL) GOTO 20
    CALL ERROR(NOR, 1)
    GOTO 8888
20      VISITL = T
  GOTO 10
!
!     *******************************
!     RUN LANCZOS. CHECK CONVERGENCE.
!     *******************************
30    IF(.NOT. (CONTIN .AND. P < POPT .AND. CNEG < REST)) GOTO 40
  IF(.NOT. LANCZO(W(ALPHA), W(BETA), TNORM, W, IW)) GOTO 8888
  IF(.NOT. TRIDIG(CONTIN, FINAL, POSNU, F, W, IW)) GOTO 8888
  GOTO 30
!
40    IF(.NOT. (CONTIN .AND. CNEG == 0)) GOTO 50
  IF(.NOT. LANCZO(W(ALPHA), W(BETA), TNORM, W, IW)) GOTO 8888
  IF(.NOT. TRIDIG(CONTIN, FINAL, POSNU, F, W, IW)) GOTO 8888
  GOTO 40
!
50    IF(.NOT. (.NOT. CONTIN .AND. CNEG == 0)) GOTO 60
!       *************************************
!       CAN NOT CONTINUE, AND NO CONVERGENCE.
!       *************************************
  CALL ERROR(NOR, 2)
  GOTO 8888
!
!     **********************
!     NO MORE LANCZOS STEPS.
!     **********************
60    FINAL = T
IF(.NOT. TRIDIG(CONTIN, FINAL, POSNU, F, W, IW)) GOTO 8888
CALL FREEID(MV + P + 1)
!
!     **************************************************
!     NEW REST. REST < 0 SHOULD NEVER HAPPEN, BUT ...
!     **************************************************
REST = REST - CNEG
IF(.NOT. REST < 0) GOTO 70
  CALL ERROR(NOR, 3)
  GOTO 8888
!
70    CPALSO = F
!     *******************************************************
!     MXNEG = START OF REST MX-VECTORS.
!     IF REST = 0, WE ARE READY, AND DO NOT HAVE TO SAVE ANY.
!     COMPUTE AND SAVE REST EIGENPAIRS. NOTE CPALSO = F.
!     UPDATE NUMR.
!     *******************************************************
MXNEG = MXRST + NUMR
IF(REST == 0) MXNEG = NIL
!
IF(.NOT. SAVEXL(W(BETAPI), W(LAMBDA), W(NU), IW(POINTR), W(S),  MXNEG, PMAX, CPALSO, W, IW)) GOTO 8888
!
NUMR = NUMR + CNEG
WRI(1) = NUMR
WRR(1) = TNORM
CALL WRINFO(NOR, 1, W, IW)
IF(REST > 0) GOTO 10
GOTO 9999
!
8888  SELDOM = F
CALL ERROR(NOR, NOR)
!
9999  CNEG = SACNEG
CPOS = SACPOS
TIME(NOR) = TIME(NOR) + SECOND() - ST
RETURN
END
