LOGICAL FUNCTION DECOMP(W, IW)
USE Numeric_Kinds_Module
!
!   $RCSfile: decomp.f,v $ $Revision: 1.3 $
!   $Date: 89/07/28 09:57:07 $
!   $State: Stable $
!
!
! **********************************************************************
!
!     PURPOSE -
!
!         TO MAKE AN ACCEPTABLE LDL(T)-DECOMPOSITION OF K-MU*M. THE
!         DECOMPOSITION CAN BE REJECTED BECAUSE,
!
!         K-MU*M OR A PRINCIPAL LEADING SUBMATRIX (VER=1) OF IT IS
!         NEARLY SINGULAR (DECIDED IN LDLT, FLAG BADMU).
!
!         THE INTERVAL (OLDMU,MU) CONTAINS TOO MANY NOT YET COMPUTED
!         EIGENVALUES, I.E. REST > MXREST.
!
!         SHOULD WE HAVE SUCH A MU AN OTHER IS CHOSEN AND WE TRY
!         AGAIN. WE MAKE ONLY ONE RETRY (THE VISIT FLAGS).
!
!     PLEASE SEE THE PROGRAMMERS GUIDE FOR INFORMATION ABOUT
!     PARAMETERS NOT EXPLAINED ABOVE, AND FOR MORE DETAILS ABOUT
!     THE FUNCTION OF THE ROUTINE.
!
!
! **********************************************************************
!
REAL(Kind=WP_Kind)  ALTMU, CK, CL, COEFF, MU, NEXTMU, OLDMU, RDUMP, SECOND, ST, ST1, TIME, TIMTQL, TLDL, &
                    TMU, TOPINV, TOPM, TPRED, TSAVE, TVECOP, W(1), A, B
INTEGER   CNEG, CONV, COPT, COUNT, CPOS, DUMMY, ERRNO, IDUMP,  ITERNO, IW(1), MXREST, N, NBADMU, NMXRST, NOR, OLCPOS,   &
                P, PFCONV, PMAX, POPT, REST, RNEW, ROLD, TCONV, TRNEW,  NUMEIG, MAXW, MAXIW
LOGICAL   BADMU, F, LDLT, T, UPDATE, USEMX, VISITD, VISITM, ZERBET, REACHB, USSMXR, USEDB
COMMON   /STLMCT/ N, ITERNO, TCONV, CNEG, CPOS, OLCPOS, RNEW, ROLD, REST, P, USEMX, ZERBET
COMMON   /STLMER/ RDUMP, ERRNO, IDUMP(2)
COMMON   /STLMIN/ A, B, NUMEIG, MAXW, MAXIW
COMMON   /STLMMU/ MU, OLDMU, NEXTMU, ALTMU
COMMON   /STLMOP/ TLDL, TOPINV, TOPM, TIMTQL, TVECOP, TPRED,  TSAVE, COEFF(4), CK, CL, CONV, PFCONV, UPDATE
COMMON   /STLMPL/ PMAX, POPT, COPT, MXREST
COMMON   /STLMST/ TIME(24), COUNT(24), NBADMU, NMXRST, DUMMY
COMMON   /STLMTS/ REACHB, USSMXR, USEDB
COMMON   /STLMTF/ T, F
!
DATA      NOR  / 3 /
!
ST = SECOND()
COUNT(NOR) = COUNT(NOR) + 1
DECOMP = T
!
!     ****************
!     SET VISIT FLAGS.
!     ****************
VISITD = F
VISITM = F
!
10    ST1 = SECOND()
COUNT(11) = COUNT(11) + 1
!     ************************
!     MAKE LDLT-DECOMPOSITION.
!     ************************
TMU = MU
IF(.NOT. LDLT(TMU, TRNEW, BADMU, W, IW)) GOTO 8888
RNEW = TRNEW
TIME(11) = TIME(11) + SECOND() - ST1
CALL WRINFO(NOR, 1, W, IW)
!     ***********************************
!     CHECK IF IT IS AN ACCEPTABLE SHIFT.
!     ***********************************
IF(.NOT. BADMU) GOTO 40
  NBADMU = NBADMU + 1
  IF(.NOT. VISITD) GOTO 20
    CALL ERROR(NOR, 1)
    GOTO 8888
!
!       ***********************
!       NO. TRY AN OTHER SHIFT.
!       ***********************
20      VISITD = T
  IF(.NOT. ITERNO == 1) GOTO 30
    IF(MU == 0.0D0) MU = -1.0D-2
    MU = MU - 1.0D-2 * ABS(MU)
    TLDL = TIME(11)
    GOTO 10
!
30      MU = MU + 1.0D-1 * (MU - OLDMU)
!       ***************************
!       CHECK IF WE HAVE REACHED B.
!       ***************************
  IF(.NOT. MU >= B) GOTO 10
    REACHB = T
    USEDB  = T
    GOTO 10
!
40    VISITD = F
!
IF(.NOT. ITERNO > 1) GOTO 9999
!       ****************
!       HOW MANY REMAIN.
!       ****************
  REST = (RNEW - ROLD) - OLCPOS
  IF(.NOT. REST > MXREST) GOTO 9999
    NMXRST = NMXRST + 1
    IF(.NOT. VISITM) GOTO 50
      CALL ERROR(NOR, 2)
      GOTO 8888
!
!         *****************************
!         TOO MANY, TRY AN OTHER SHIFT.
!         *****************************
50        MU = ALTMU
    VISITM = T
!         *****************************
!         CHECK IF WE BACK AWAY FROM B.
!         *****************************
    IF(.NOT. MU < B) GOTO 10
      REACHB = F
      USEDB  = F
      GOTO 10
!
8888  DECOMP = F
CALL ERROR(NOR, NOR)
!
9999  TIME(NOR) = TIME(NOR) + SECOND() - ST
RETURN
END
