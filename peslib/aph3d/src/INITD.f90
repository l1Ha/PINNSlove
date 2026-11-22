SUBROUTINE INITD(N, A, B, MAXL, PROFIL, PMAX, MXREST, MSGLVL, MAXW, MAXIW, &
       DAFILE, MAXREC, KFILE, X, BG, TCONV, NLEFT, ERRNO, W, IW, LEN)
USE Numeric_Kinds_Module
!
!   $RCSfile: initd.f,v $ $Revision: 1.3 $
!   $Date: 89/07/28 09:57:14 $
!   $State: Stable $
!
!
! **********************************************************************
!
!     PURPOSE -
!
!         INITALIZES ALL THE COMMON BLOCKS. SOME INITIALIZATIONS ARE
!         MADE ONLY SO THAT A PREMATURE TERMINATION WILL NOT CAUSE A
!         SUBSEQUENT ABORT WHEN TRYING TO WRITE THE VALUE OF AN
!         A NOT INITIALIZED VARIABLE. OTHERS ARE MADE BECAUSE WRINF4
!         WOULD OTHERWISE SOMETIMES WRITE VALUES NOT INITIALIZED.
!
!
!     OUTPUT PARAMETERS -
!
!         LEN HAS THE VALUE COMPUTED IN CHECK.
!
!
!         THE REASON FOR THE VARIABLE NAMES BELOW (E.G. N1) IS THAT
!         IT IS FORBIDDEN TO HAVE DUMMY PARAMETERS IN COMMON.
!
!     PLEASE SEE THE USER GUIDE FOR INFORMATION ABOUT
!     PARAMETERS NOT EXPLAINED ABOVE.
!
!
! **********************************************************************
!
REAL(Kind=WP_Kind)  A, A1, ALTMU, B, B1, CK, CL, COEFF, FACTOR, SRELPR, MU,  NEXTMU, OLDMU, RDUMP, SMALL, TEMP, TIME, TIMTQL, &
                TLDL, TOLBPI, TOLLDL, TOLPDM, TOLS1I, TOLZBT, TOLZNU, TOPINV, TOPM, TPRED, TSAVE, TVECOP, W(1), WRR, BG
INTEGER   ACTIVE, ADDRSS, ALPHA, BETA, BETAPI, CNEG, CNEGF, CONV,  COPT, COUNT, CPOS, D, DAFIL1, DAFILE, DUMMY, ERRNO,         &
                FREE, I, I1, IDUMP, ITERNO, IW(1), K, KFILE, KFILE1, LAMBDA, LEFTP, LEN, LENADR, LP, M, MAXIW, MAXIW1, MAXL,   &
                MAXL1, MAXRE1, MAXREC, MAXW, MAXW1, MV, MXNEW, MXOLD, MXREST, MXRST, N, N1, N2, NBADMU, NIL, NMXRST,           &
                NOACTN, NREAD, NU, NUMEIG, NUMVEC, NWRITE, OLCPOS, P, PFCONV, PMAX, POINTR, POPT, PV1, PV2,                    &
                READID, READK, REST, RFIRST, RIGHTC, RIGHTM, RIGHTP,  RNEW, ROLD, S, SAEVAL, SAVE, SAVFRE, SCPX, SCR, SOLCPX,  &
                STADEW, TCONV, V, VER, VTEMP, WAD, WRI, WRITID, X, TCONV1, NLEFT, ERRNO1, MSGLV1, MSGLVL, PMAX1,               &
                PROFIL, PROFI1, MXRES1, NERR  , NOUT , X1
LOGICAL   DIAGM, F, MEQI, SAFRST, T, UPDATE, USEMX, WRL, ZERBET, REACHB, USSMXR, USEDB
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
COMMON   /STLMVR/ SMALL, VER
COMMON   /STLMWH/ WAD(2), ACTIVE(2)
COMMON   /STLMWR/ WRR(5), WRI(5), WRL(5)
!
!     *****************************
!     INITIALIZE OUTPUT PARAMETERS.
!     *****************************
X = 0
BG = A
WRR(1) = A
TCONV = 0
NLEFT = 0
ERRNO = 0
WRI(1) = 0
!     ************
!     CHECK INPUT.
!     ************
COUNT(4) = 0
CALL CHECKINP(N, A, B, MAXL, PROFIL, PMAX, MXREST, MSGLVL, MAXW, MAXIW, DAFILE, MAXREC, KFILE, ERRNO, W, IW, LEN)
!      Write(Msg_Unit,*)'AFTER CHECK IN INITD,LEN=', LEN
IF(LEN > 0) GOTO 9999
!
!
!     *************************
!     INITIALIZE COMMON BLOCKS.
!     *************************
!
!
!     **********
!     / STLMAC /
!     **********
NOACTN = 0
FREE = 1
SAVE = 2
SAVFRE = 3
!     **********
!     / STLMCT /
!     **********
N1 = N
ITERNO = 0
TCONV1 = 0
CNEG = 0
CPOS = 0
OLCPOS = 0
RNEW = 0
ROLD = 0
REST = 0
P = 0
USEMX = .TRUE.
ZERBET = .FALSE.
!     **********
!     / STLMDS /
!     **********
K = 1
M = 1
D = 1
DIAGM = PROFIL >= 2
!     **********
!     / STLMMI /
!     **********
MEQI = PROFIL == 3
!     ****************************************
!     PV1 POINTS TO THE START OF ALPHA IN W.
!     PV2 POINTS TO THE START OF POINTR IN IW.
!     ****************************************
PV1 = 1
PV2 = 1
IF(.NOT. VER == 1) GOTO 10
  M = IW(N) + 1
  IF(MEQI) PV1 = M
  IF(DIAGM .AND. .NOT. MEQI) PV1 = M + N
  IF(.NOT. DIAGM) PV1 = M + IW(N)
  PV2 = N + 1
!       *******
!       SAVE K.
!       *******
  REWIND KFILE
  I1 = IW(N)
  WRITE(KFILE) (W(I), I=1, I1)
!     **********
!     / STLMER /
!     **********
10    RDUMP = 0.0D0
ERRNO1 = 0
IDUMP(1) = 0
IDUMP(2) = 0
!     **********
!     / STLMFT /
!     **********
CNEGF = 0
RFIRST = 0
SAFRST = .FALSE.
!     **********
!     / STLMID /
!     **********
NIL = 0
MV = 0
V = 0
MXNEW = 0
MXOLD = 0
MXRST = 0
SCPX = 0
SOLCPX = 0
X1 = 0
!     **********
!     / STLMIO /
!     **********
SAEVAL = 1
READID = 2
WRITID = 3
READK = 4
N2 = N
!     **********
!     / STLMMU /
!     **********
MU = A
OLDMU = 0.0D0
NEXTMU = 0.0D0
ALTMU = 0.0D0
!     **********
!     / STLMOP /
!     **********
TOPINV =  -1.0D0
TOPM =  -1.0D0
TIMTQL =  -1.0D0
TVECOP =  -1.0D0
TPRED = 0.0D0
TSAVE = 0.0D0
DO 20 I = 1, 4
  COEFF(I) = 0.0D0
20      CONTINUE
CK = 2.0D0
CL = 1.0D1
CONV = 0
PFCONV = 0
UPDATE = .TRUE.
TLDL = 0.0D0
!     ***********
!     / STLMPL /
!     **********
POPT = 0
COPT = 0
MXRES1 = MXREST
PMAX1 = PMAX
!     **************************************
!     / STLMPR / SEE BLOCKDATA PROGRAM ALSO.
!     **************************************
MSGLV1 = MSGLVL
!     **********
!     / STLMPF /
!     **********
PROFI1 = PROFIL
!     **********
!     / STLMPV /
!     **********
ALPHA = PV1
BETA = ALPHA + PMAX
BETAPI = BETA + PMAX
LAMBDA = BETAPI + PMAX
NU = LAMBDA + PMAX
POINTR = PV2
S = NU + PMAX
SCR = BETAPI
!     **********
!     / STLMST /
!     **********
DO 30 I = 1, 24
  TIME(I) = 0.0D0
  COUNT(I) = 0
30      CONTINUE
NBADMU = 0
NMXRST = 0
!     **********
!     / STLMTF /
!     **********
T = .TRUE.
F = .FALSE.
!     **************************************
!     / STLMTL / SEE BLOCKDATA PROGRAM ALSO.
!     **************************************
!
TOLBPI = 1.0D-7
TEMP = N
TOLS1I = 1.0D-2 / SQRT(TEMP)
FACTOR = 1.0D-6
TOLPDM = SMALL * TEMP
TOLZBT = 1.0D-5
TOLZNU = SMALL * 1.0D5
TOLLDL(1) = SMALL * 1.0D5
TOLLDL(2) = 1.0D5
TOLLDL(3) = -1.0D0
!     **********
!     / STLMTS /
!     **********
REACHB = .FALSE.
USSMXR = MXREST >= 1
USEDB  = .FALSE.
!     **********
!     / STLMUI /
!     **********
!     **********
!     / STLMWH /
!     **********
ACTIVE(1) = NIL
ACTIVE(2) = NIL
WAD(1) = S + PMAX**2
WAD(2) = WAD(1) + N
!     **********
!     / STLMWR /
!     **********
DO 50 I = 1, 5
  WRR(I) = 0.0D0
  WRI(I) = 0
  WRL(I) = .TRUE.
50      CONTINUE
!
!
!     **********
!     / STLMAD /
!     **********
ADDRSS = POINTR + PMAX + 1
IF(VER == 3) ADDRSS = 0
DAFIL1 = DAFILE
IF(VER == 3) DAFIL1 = 0
KFILE1 = KFILE
IF(VER >= 2) KFILE1 = 0
LP = WAD(2) + N
IF(VER == 3) LP = WAD(1)
MAXL1 = MAXL
NREAD = 0
NWRITE = 0
!
!     **********
!     / STLMEW /
!     ***********
STADEW = LP + MAXL
IF(.NOT. VER == 3) GOTO 60
  STADEW = 0
  NUMVEC = 0
  LENADR = 0
  MAXRE1 = 0
  GOTO 80
!
!
!
60    NUMVEC = (MAXW - STADEW + 1) / N
IF(NUMVEC < 0) NUMVEC = 0
LENADR = 5 * PMAX + 2 + MAXL
IF(MEQI) LENADR = PMAX + 1 + MAXL
!
I1 = ADDRSS + LENADR - 1
DO 70 I = ADDRSS, I1
  IW(I) = 0
70      CONTINUE
!
!     *****************************************************
!     SINCE USEMX ETC. ARE NOT KNOWN AT THIS STAGE, WE ONLY
!     GIVE ADDRESSES TO V + 1 (AND MV + 1), SINCE THEY MAY
!     BE USED IN FRSTIT.
!     *****************************************************
VTEMP = 0
IF(.NOT. MEQI) VTEMP = PMAX + 1
I1 = ADDRSS + VTEMP
IF(NUMVEC >= 1) IW(I1) = -STADEW
I1 = ADDRSS + MV
IF(NUMVEC >= 2 .AND. .NOT. MEQI) IW(I1) = -STADEW - N
!
LEFTP = 1
RIGHTP = MAXREC
RIGHTC = MAXREC + 1
RIGHTM = RIGHTC
MAXRE1 = MAXREC
!
!
!     **********
!     / STLMIN /
!     **********
80    A1 = A
B1 = B
NUMEIG = 0
MAXW1 = MAXW
MAXIW1 = MAXIW
!
!     *********************************
!     / STLMVR / SEE BLOCKDATA PROGRAM.
!     *********************************
!
!
9999  RETURN
END
