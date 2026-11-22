LOGICAL FUNCTION LANCZO(ALPHA, BETA, TNORM, W, IW)
USE Numeric_Kinds_Module
!
!   $RCSfile: lanczo.f,v $ $Revision: 1.3 $
!   $Date: 89/07/28 09:57:18 $
!   $State: Stable $
!
!
! **********************************************************************
!
!     PURPOSE -
!
!         TAKES ONE LANCZOS STEP.
!
!
!     INPUT PARAMETERS -
!
!         TNORM = MAXIMUM NORM OF THE TRIDIAGONAL T-MATRIX OF
!         DIMENSION P (> 1).
!
!
!     OUTPUT PARAMETERS -
!
!         TNORM = UPDATED TNORM, I.E. NORM OF T OF ORDER P+1.
!
!     PLEASE SEE THE PROGRAMMERS GUIDE FOR INFORMATION ABOUT
!     PARAMETERS NOT EXPLAINED ABOVE, AND FOR MORE DETAILS ABOUT
!     THE FUNCTION OF THE ROUTINE.
!
!
! **********************************************************************
!
REAL(Kind=WP_Kind)  ALPHA(1), BETA(1), DOT, FACTOR, SRELPR, RDUMP, SECOND, ST, ST1, TIME, TNORM, TOLBPI, TOLLDL, TOLPDM, TOLS1I,   &
                TOLZBT, TOLZNU, W(1), XX
INTEGER   ACTION, CNEG, COUNT, CPOS, DUMMY, ERRNO, FREE, I, IDUMP, ITERNO, IW(1), MV, MXNEW, MXOLD, MXRST, N, NBADMU,            &
                NEXTV, NIL, NMXRST, NOACTN, NOR, OLCPOS, P, REST, RNEW, ROLD, SAVE, SAVFRE, SCPX, SOLCPX, TCONV, V, X
LOGICAL   F, MEQI, MULVEC, OPINV, OPM, SCPROD, SUBVEC, T, USEMX, ZERBET
COMMON   /STLMAC/ NOACTN, FREE, SAVE, SAVFRE
COMMON   /STLMCT/ N, ITERNO, TCONV, CNEG, CPOS, OLCPOS, RNEW, ROLD, REST, P, USEMX, ZERBET
COMMON   /STLMER/ RDUMP, ERRNO, IDUMP(2)
COMMON   /STLMID/ NIL, MV, V, MXNEW, MXOLD, MXRST, SCPX, SOLCPX, X
COMMON   /STLMMI/ MEQI
COMMON   /STLMST/ TIME(24), COUNT(24), NBADMU, NMXRST, DUMMY
COMMON   /STLMTF/ T, F
COMMON   /STLMTL/ SRELPR, TOLBPI, TOLS1I, FACTOR, TOLPDM, TOLZBT, TOLZNU, TOLLDL(3)
!
DATA      NOR  / 10 /
!
ST = SECOND()
COUNT(NOR) = COUNT(NOR) + 1
LANCZO = T
!
P = P + 1
!     ********************************
!     NEXTV = ID TO THE NEXT V-VECTOR.
!     ********************************
NEXTV = V + P + 1
ACTION = FREE
IF(MEQI) ACTION = NOACTN
!
IF(.NOT. OPINV(NEXTV, MV + P, W, IW)) GOTO 8888
!
IF(.NOT. P > 1) GOTO 10
  CALL FREEID(MV + P)
  IF(.NOT. SUBVEC(NEXTV, V + P - 1, BETA(P - 1), W, IW, FREE, 1))GOTO 8888
!
10    IF(.NOT. SCPROD(NEXTV, MV + P, ALPHA(P), W, IW, ACTION, 1))GOTO 8888
!
IF(.NOT. SUBVEC(NEXTV, V + P, ALPHA(P), W, IW, FREE, 1)) GOTO 8888
!
!     ********************
!     REORTHOGONALIZATION.
!     ********************
ST1 = SECOND()
DO 20 I = 1, P
  COUNT(24) = COUNT(24) + 1
  IF(.NOT. SCPROD(NEXTV, MV + I, DOT, W, IW, ACTION, 2)) GOTO 8888
  IF(.NOT. SUBVEC(NEXTV, V + I, DOT, W, IW, FREE, 3)) GOTO 8888
20      CONTINUE
TIME(24) = TIME(24) + SECOND() - ST1
!
IF(MEQI) GOTO 25
IF(.NOT. OPM(MV + P + 1, NEXTV, W, IW)) GOTO 8888
25    IF(.NOT. SCPROD(MV + P + 1, NEXTV, DOT, W, IW, NOACTN, 1))GOTO 8888
!
!
IF(.NOT. P > 2) GOTO 30
!       *************
!       UPDATE TNORM.
!       *************
  TNORM = MAX(TNORM, BETA(P - 2) + ABS(ALPHA(P - 1)) +BETA(P - 1))
  GOTO 40
!
30    IF(P == 1) TNORM = ABS(ALPHA(1))
IF(P == 2) TNORM = TNORM + BETA(1)
!
40    BETA(P) =  SIGN(SQRT(ABS(DOT)), DOT)
!
IF(.NOT. ABS(BETA(P)) <= TNORM * TOLZBT) GOTO 50
!       ***********************************************
!       BETA(P) IS CONSIDERED TO BE ZERO. HALT LANCZOS.
!       ***********************************************
  XX = 1.0D0
  ZERBET = T
  GOTO 60
!
50    IF(.NOT. BETA(P) <  (- TNORM * TOLZBT) ) GOTO 60
!       *********************************
!       M IS CONSIDERED TO BE INDEFINITE.
!       *********************************
  RDUMP = DOT
  CALL ERROR(NOR, 1)
  GOTO 8888
!
60    IF(.NOT. ZERBET) XX = 1.0D0 / BETA(P)
ACTION = SAVFRE
IF(MEQI) ACTION = SAVE
IF(.NOT. MULVEC(NEXTV, NEXTV, XX, W, IW, ACTION, 1)) GOTO 8888
IF(MEQI) GOTO 9999
IF(.NOT. MULVEC(MV + P + 1, MV + P + 1, XX, W, IW, SAVE, 1))GOTO 8888
!
GOTO 9999
!
8888  LANCZO = F
CALL ERROR(NOR, NOR)
!
9999  TIME(NOR) = TIME(NOR) + SECOND() - ST
RETURN
END
