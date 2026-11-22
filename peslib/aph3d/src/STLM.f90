SUBROUTINE STLM(N, A, B, MAXL, PROFIL, PMAX, MXREST, MSGLVL, MAXW, MAXIW, DAFILE, MAXREC, KFILE, X, BG, TCONV, NLEFT, ERRNO, W, IW)
USE Numeric_Kinds_Module
!
!   $RCSfile: stlm.f,v $ $Revision: 1.3 $
!   $Date: 89/07/28 09:57:44 $
!   $State: Stable $
!
!
! **********************************************************************
!
!     PURPOSE -
!
!         MAIN ROUTINE FOR STLM.
!
!         PLEASE SEE THE USER GUIDE FOR FURTHER COMMENTS.
!
!
! **********************************************************************
!
!
SAVE
REAL(Kind=WP_Kind)  A, A1, ALTMU, B, B1, CK, CL, COEFF, FACTOR, SRELPR, MU,  NEXTMU, OLDMU, RDUMP, RV, SECOND, SMALL, TIME, TIMTQL
REAL(Kind=WP_Kind)  TLDL, TOLBPI, TOLLDL, TOLPDM, TOLS1I, TOLZBT, TOLZNU, TOPINV, TOPM, TPRED, TSAVE, TVECOP, W(1), WRR, BG
INTEGER   ACTIVE, ADDRSS, ALPHA, BETA, BETAPI, CNEG, CNEGF, CONV, COPT, COUNT, CPOS, D, DAFIL1, DAFILE, DUMMY, ERRNO
INTEGER   FREE, IDUMP, ITERNO, IV, IW(1), K, KFILE, KFILE1,LAMBDA, LEFTP, LEN, LENADR, LP, M, MAXIW, MAXIW1, MAXL
INTEGER   MAXL1, MAXRE1, MAXREC, MAXW, MAXW1, MSGLVL, MV, MXNEW, MXOLD, MXREST, MXRST, N, N1, N2, NBADMU, NIL, NMXRST
INTEGER   NOACTN, NOR, NREAD, NU, NUMEIG, NUMVEC, NWRITE,  OLCPOS, P, PFCONV, PMAX, POINTR, POPT, READID, READK
INTEGER   REST, RFIRST, RIGHTC, RIGHTM, RIGHTP, RNEW, ROLD, S, SAEVAL, SAVE, SAVFRE, SCPX, SCR, SOLCPX, STADEW, TCONV
INTEGER   NERR  , NOUT , V, VER, WAD, WRI, WRITID, X, X1,  TCONV1, NLEFT, ERRNO1, PMAX1, MXRES1, MSGLV1, PROFIL
INTEGER   PROFI1
LOGICAL   ALWAYS, DECOMP, DIAGM, F, FRSTIT, MEQI,  SAFRST, SELDOM, T, TERMIN, UPDATE, USEMX, WRL, ZERBET,  REACHB, USSMXR, USEDB
COMMON   /STLMAC/ NOACTN, FREE, SAVE, SAVFRE
COMMON   /STLMAD/ ADDRSS, DAFIL1, KFILE1, LP, MAXL1, NREAD, NWRITE
COMMON   /STLMCT/ N1, ITERNO, TCONV1, CNEG, CPOS, OLCPOS, RNEW, ROLD, REST, P, USEMX, ZERBET
COMMON   /STLMDS/ K, M, D, DIAGM
COMMON   /STLMER/ RDUMP, ERRNO1, IDUMP(2)
COMMON   /STLMEW/ LEFTP, LENADR, MAXRE1, NUMVEC, RIGHTC, RIGHTM, RIGHTP, STADEW
COMMON   /STLMFT/ CNEGF, RFIRST, SAFRST
COMMON   /STLMID/ NIL, MV, V, MXNEW, MXOLD, MXRST, SCPX, SOLCPX, X1
COMMON   /STLMIN/ A1, B1, NUMEIG, MAXW1, MAXIW1
COMMON   /STLMIO/ SAEVAL, READID, WRITID, READK, N2
COMMON   /STLMMI/ MEQI
COMMON   /STLMMU/ MU, OLDMU, NEXTMU, ALTMU
COMMON   /STLMOP/ TLDL, TOPINV, TOPM, TIMTQL, TVECOP, TPRED, TSAVE, COEFF(4), CK, CL, CONV, PFCONV, UPDATE
COMMON   /STLMPL/ PMAX1, POPT, COPT, MXRES1
COMMON   /STLMPF/ PROFI1
COMMON   /STLMPR/ MSGLV1, NERR  , NOUT
COMMON   /STLMPV/ ALPHA, BETA, BETAPI, LAMBDA, NU, POINTR, S, SCR
COMMON   /STLMST/ TIME(24), COUNT(24), NBADMU, NMXRST, DUMMY
COMMON   /STLMTF/ T, F
COMMON   /STLMTL/ SRELPR, TOLBPI, TOLS1I, FACTOR, TOLPDM, TOLZBT, TOLZNU, TOLLDL(3)
COMMON   /STLMTS/ REACHB, USSMXR, USEDB
COMMON   /STLMUI/ RV(5), IV(9)
COMMON   /STLMVR/ SMALL, VER
COMMON   /STLMWH/ WAD(2), ACTIVE(2)
COMMON   /STLMWR/ WRR(5), WRI(5), WRL(5)
!
DATA      NOR  / 20 /
!
!     ********************
!     INITIALIZATION PART.
!     ********************
CALL INITD(N, A, B, MAXL, PROFIL, PMAX, MXREST, MSGLVL, MAXW, MAXIW, DAFILE, MAXREC, KFILE, X, BG, TCONV, NLEFT, ERRNO, W, IW, LEN)
IF(LEN > 0) RETURN
!
!     ***********************
!     WRITE INPUT PARAMETERS.
!     ***********************
CALL WRINFO(NOR, 1, W, IW)
!
TIME(NOR) = SECOND()
COUNT(NOR) = COUNT(NOR) + 1
!
!     **********
!     MAIN LOOP.
!     **********
10      ITERNO = ITERNO + 1
  TSAVE = SECOND()
  CALL WRINFO(NOR, 2, W, IW)
!       *******************************************
!       COMPUTE LDL(T)-DECOMPOSITION OF K - MU * M.
!       *******************************************
  IF(.NOT. DECOMP(W, IW)) GOTO 8888
!
  CALL WRINFO(NOR, 3, W, IW)
  IF(USEDB .AND. REST == 0) GOTO 6666
!
  IF(.NOT. ITERNO == 1) GOTO 20
!         ****************************************
!         INITIALIZATIONS FOR THE FIRST ITERATION.
!         ****************************************
    RV(1) = MU
    IF(.NOT. FRSTIT(W, IW)) GOTO 8888
!
!       *************************************
!       HANDLES THE COMPUTATIONS FOR A SHIFT.
!       *************************************
20      IF(VER < 3) CALL INITEW(IW(ADDRSS))
  IF(.NOT. ALWAYS(W, IW)) GOTO 8888
!
  IF(.NOT. REST > 0) GOTO 30
!         ********************
!         MISSING EIGENVALUES.
!         ********************
    IF(.NOT. SELDOM(W, IW)) GOTO 8888
!
30      CALL WRINFO(NOR, 4, W, IW)
!       *********************************************
!       CHECK IF WE HAVE REACHED THE LAST EIGENVALUE.
!       *********************************************
  IF(CPOS + RNEW == N) GOTO 7777
!
!       ************************************
!       DECIDE IF WE SHOULD CONTINUE OR NOT.
!       ************************************
6666    IF(TERMIN(1, W, IW, BG, TCONV, NLEFT, ERRNO)) GOTO 9999
!       *********************************
!       UPDATE CERTAIN CONTROL VARIABLES.
!       *********************************
  CALL REPL
!
  GOTO 10
!
!     ************************************
!     WE HAVE REACHED THE LAST EIGENVALUE.
!     ************************************
7777  IF(TERMIN(2, W, IW, BG, TCONV, NLEFT, ERRNO)) GOTO 9999
GOTO 9999
!
!     ***********************************
!     ERROR EXIT, UNLESS ERRORS IN INPUT.
!     ***********************************
8888  CALL ERROR(NOR, NOR)
IF(TERMIN(3, W, IW, BG, TCONV, NLEFT, ERRNO)) GOTO 9999
!
9999  TIME(NOR) = SECOND() - TIME(NOR)
!
!
!     *************************
!     WRITE OUTPUT INFORMATION.
!     *************************
CALL WRINFO(NOR, 5, W, IW)
!
X = X1
!
RETURN
END
