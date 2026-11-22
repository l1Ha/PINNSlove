LOGICAL FUNCTION LDLT(MU, RNEW, BADMU, W, IW)
USE Numeric_Kinds_Module
!
!   $RCSfile: ldlt.f,v $ $Revision: 1.3 $
!   $Date: 89/07/28 09:57:22 $
!   $State: Stable $
!
!
! **********************************************************************
!
!     PURPOSE -
!
!         DRIVER ROUTINE FOR LDL(T)-FACTORIZATION ROUTINES.
!
!
!     OUTPUT PARAMETERS -
!
!         BADMU = TRUE, IF MU IS CONSIDERED TO BE UNSUITABLE.
!
!     PLEASE SEE THE PROGRAMMERS GUIDE FOR INFORMATION ABOUT
!     PARAMETERS NOT EXPLAINED ABOVE, AND FOR MORE DETAILS ABOUT
!     THE FUNCTION OF THE ROUTINE.
!
!
! **********************************************************************
!
REAL(Kind=WP_Kind) MU, SMALL, W(1)
INTEGER   ACTIVE, AD, ADDRSS, D, DAFILE, IW(1), K, KFILE, LEN, LP, M, MAXL, N, NOR, NREAD, NUMEL, NWRITE, READID, READK,    &
                RNEW, SAEVAL, VER, WAD1, WAD2, WRITID, PROFIL
LOGICAL   BADMU, DIAGM, F, IO, MEQI, T
COMMON   /STLMAD/ ADDRSS, DAFILE, KFILE, LP, MAXL, NREAD, NWRITE
COMMON   /STLMDS/ K, M, D, DIAGM
COMMON   /STLMIO/ SAEVAL, READID, WRITID, READK, N
COMMON   /STLMMI/ MEQI
COMMON   /STLMPF/ PROFIL
COMMON   /STLMTF/ T, F
COMMON   /STLMVR/ SMALL, VER
COMMON   /STLMWH/ WAD1, WAD2, ACTIVE(2)
!
DATA      NOR  / 11 /
!
LDLT = T
!
IF(.NOT. VER == 1) GOTO 10
  AD = D + N - 1
!       *******************************
!       FETCH K FROM SECONDARY STORAGE.
!       *******************************
  NUMEL = IW(AD)
  IF(.NOT. IO(W(K), 1, READK, NUMEL, IW(ADDRSS) )) GOTO 8888
!
!       *****************************************************
!       MAKE DECOMPOSITION. WAD1 AND WAD2 ARE ADDRESSES TO
!       N-VECTORS IN W. THEY ARE USED AS SCRATCH VECTORS AND
!       DO NOT HAVE ANY IDENTIFIERS. (WE DO NOT NEED ALLOC IN
!       THIS SPECIAL CASE).
!       *****************************************************
  CALL LDLSUB(W(K), W(M), IW(D), W(WAD1), W(WAD2), MU, N, RNEW, DIAGM, MEQI, LEN)
!
  BADMU = F
  IF(LEN == (-1)) BADMU = T
  GOTO 9999
!
10    IF(VER == 3) CALL LDL3(W, IW, MU, N, RNEW, PROFIL, LEN)
!
IF(VER == 2) CALL LDL2(W, IW, W(WAD1), W(WAD2), MU, N, RNEW, PROFIL, LEN)
!
BADMU = F
IF(LEN == (-1)) BADMU = T
IF(.NOT. LEN > 0) GOTO 9999
  CALL ERROR(100, LEN)
  CALL ERROR(100, 100)
!
8888  LDLT = F
CALL ERROR(NOR, NOR)
!
9999  RETURN
END
