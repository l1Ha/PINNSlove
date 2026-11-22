SUBROUTINE CHECKINP(N, A, B, MAXL, PROFIL, PMAX, MXREST, MSGLVL, MAXW, MAXIW, DAFILE, MAXREC, KFILE,  ERRNO, W, IW, LEN)
USE Numeric_Kinds_Module
SAVE
!
!   $RCSfile: check.f,v $ $Revision: 1.3 $
!   $Date: 89/07/28 09:57:04 $
!   $State: Stable $
!
!
! **********************************************************************
!
!     PURPOSE -
!
!         THIS ROUTINE CHECKS INPUT.
!
!
!     OUTPUT PARAMETERS -
!
!         LEN =  0, NO ERROR HAS BEEN DETECTED.
!                1, ILLEGAL VER.
!                2, N LE 0.
!                3, ILLEGAL PROFIL.
!                4, MAXL LE 0.
!                5, MAXREC LE 0.
!                6, ERROR IN D (VER = 1 ONLY).
!                7, ERROR IN M (VER = 1 ONLY).
!                8, ILLEGAL MSGLVL.
!                9, ILLEGAL PMAX.
!               10, TOO SMALL MAXW.
!               11, TOO SMALL MAXIW.
!               12, B LE A.
!               13, ILLEGAL FILE UNIT(S).
!
!     PLEASE SEE THE USER GUIDE FOR INFORMATION ABOUT
!     PARAMETERS NOT EXPLAINED ABOVE, AND FOR MORE DETAILS ABOUT
!     THE FUNCTION OF THE ROUTINE.
!
!
! **********************************************************************
!
REAL(Kind=WP_Kind) A, A1, B, B1, MU, RDUMP, SMALL, TOTALT, TPRED, W(1)
INTEGER   ADDRSS, CNEGF, COPT, CPOS, DAFIL1, DAFILE, ERRNO,  ERRNO1, I, IDUM, IDUMP, ITERNO, IW(1), KFILE, KFILE1,     &
          LEFTP, LEN, LENADR, LP, MAXIW, MAXIW1, MAXIWT, MAXL, MAXL1, MAXRE1, MAXREC, MAXW, MAXW1, MAXWT, MDIM,        &
          MSGLVL, MXREST, N, N1, N2, NREAD, NUMEIG,  NUMVEC, NWRITE, PMAX, POPT, RFIRST, RIGHTC, RIGHTM,               &
          RIGHTP, RNEW, STADEW, TCONV, NERR  , NOUT , VER, X,   ERRNO2, PROFI1, PROFIL, MXRES1, PMAX1, MSGLV1
LOGICAL   LDUM, DIAGM, MEQI
COMMON   /STLMAD/ ADDRSS, DAFIL1, KFILE1, LP, MAXL1, NREAD, NWRITE
COMMON   /STLMCT/ N2, IDUM(9), LDUM(2)
COMMON   /STLMER/ RDUMP, ERRNO1, IDUMP(2)
COMMON   /STLMEW/ LEFTP, LENADR, MAXRE1, NUMVEC, RIGHTC, RIGHTM, RIGHTP, STADEW
COMMON   /STLMIN/ A1, B1, NUMEIG, MAXW1, MAXIW1
COMMON   /STLMPF/ PROFI1
COMMON   /STLMPL/ PMAX1, POPT, COPT, MXRES1
COMMON   /STLMPR/ MSGLV1, NERR  , NOUT
COMMON   /STLMUI/ MU(3), TOTALT, TPRED, CNEGF, CPOS, ERRNO2, ITERNO, N1, RFIRST, RNEW, TCONV, X
COMMON   /STLMVR/ SMALL, VER
!
!
DIAGM = PROFIL >= 2
MEQI  = PROFIL == 3
!
LEN = 1
IF(VER < 1 .OR. VER > 3) GOTO 8888
!
LEN = 2
IF(N <= 0) GOTO 8888
!
LEN = 3
IF(PROFIL < 1  .OR.  PROFIL > 3) GOTO 8888
!
LEN = 4
IF(MAXL <= 0) GOTO 8888
!
IF(.NOT. VER <= 2) GOTO 5
  LEN = 5
  IF(MAXREC <= 0) GOTO 8888
!
5     IF(.NOT. VER == 1) GOTO 10
  CALL DCHECK(IW, N, LEN, IDUMP(1))
  IF(LEN == 6) GOTO 8888
!
  I = IW(N) + 1
  CALL MCHECK(W(I), IW, N, DIAGM, MEQI, LEN, IDUMP(1))
  IF(LEN == 7) GOTO 8888
!
10    LEN = 8
IF(MSGLVL < 0 .OR. MSGLVL > 4) GOTO 8888
!
LEN = 9
IF(PMAX < 15) GOTO 8888
!
IF(.NOT. VER == 1) GOTO 20
  MDIM = IW(N)
  IF(DIAGM) MDIM = N
  IF(MEQI)  MDIM = 0
  MAXWT = IW(N) + MDIM + PMAX * (5 + PMAX) + 2 * N + MAXL
  MAXIWT = N + 6 * PMAX + MAXL + 3
  IF(MEQI) MAXIWT = N + 2 * PMAX + MAXL + 2
  GOTO 40
!
20    IF(.NOT. VER == 2) GOTO 30
  MAXWT = PMAX * (5 + PMAX) + 2 * N + MAXL
  MAXIWT = 6 * PMAX + MAXL + 3
  IF(MEQI) MAXIWT = 2 * PMAX + MAXL + 2
  GOTO 40
!
30    MAXWT = PMAX * (5 + PMAX) + MAXL
MAXIWT = PMAX + 1 + MAXL
!
40    LEN = 10
IF(MAXW < MAXWT)THEN
Write(Msg_Unit,*) 'maxw<maxwt in check', maxw, maxwt
Write(Msg_Unit,*) 'pmax, n, maxl, mdim=', pmax, n, maxl, mdim

GOTO 8888
ENDIF
!
LEN = 11
IF(MAXIW < MAXIWT) GOTO 8888
!
LEN = 12
IF(B <= A) GOTO 8888
!
LEN = 13
IF(.NOT. VER <= 2 ) GOTO 8000
  IF(DAFILE == NERR    .OR.  DAFILE == NOUT ) GOTO 8888
  IF(VER == 1  .AND.  (DAFILE == KFILE .OR. KFILE == NERR    .OR.  KFILE == NOUT ) ) GOTO 8888
!
!
8000  LEN = 0
GOTO 9999
!
!     *************************************************
!     ERROR EXIT.
!     FILL STLMUI AND THE COMMON BLOCKS USED IN WRINFO.
!     *************************************************
8888  MU(1) = A
MU(2) = A
MU(3) = A
TOTALT = 0.0D0
TPRED = 0.0D0
CNEGF = 0
CPOS = 0
ITERNO = 0
N1 = N
RFIRST = 0
RNEW = 0
TCONV = 0
X = 0
!
N2 = N
A1 = A
B1 = B
MAXW1 = MAXW
MAXIW1 = MAXIW
MAXL1 = MAXL
KFILE1 = KFILE
DAFIL1 = DAFILE
MAXRE1 = MAXREC
PROFI1 = PROFIL
PMAX1 = PMAX
MXRES1 = MXREST
MSGLV1 = 0
IF(MSGLVL == 1) MSGLV1 = 1
IF(MSGLVL > 1 .AND. MSGLVL <= 4) MSGLV1 = 2
CALL WRINFO(20, 1, W, IW)
!
CALL ERROR(25, LEN)
CALL ERROR(25, 25)
CALL ERROR(20, 20)
ERRNO = ERRNO1
ERRNO2 = ERRNO1
!
!
9999  RETURN
END
