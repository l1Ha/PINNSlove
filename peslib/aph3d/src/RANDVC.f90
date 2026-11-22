LOGICAL FUNCTION RANDVC(ID, W, IW, SAVE1)
USE Numeric_Kinds_Module
!
!   $RCSfile: randvc.f,v $ $Revision: 1.3 $
!   $Date: 89/07/28 09:57:34 $
!   $State: Stable $
!
!
! **********************************************************************
!
!     PURPOSE -
!
!         VEC(ID) = (K - MU * M)**(-1) * M * RANDOMVECTOR.
!
!
!     INPUT PARAMETERS -
!
!         ID    = IDENTIFIER OF THE RESULTING VECTOR.
!         SAVE1 = TRUE, VEC(ID) IS STORED AFTER INITIALIZATION.
!                 FALSE, VEC(ID) IS NOT STORED.
!                 SAVE1 IS ONLY USED IF VER = 1 OR 2.
!
!     PLEASE SEE THE PROGRAMMERS GUIDE FOR INFORMATION ABOUT
!     PARAMETERS NOT EXPLAINED ABOVE, AND FOR MORE DETAILS ABOUT
!     THE FUNCTION OF THE ROUTINE.
!
!
! **********************************************************************
!
REAL(Kind=WP_Kind)  RANFX, SECOND, ST, TIME, W(1), SMALL
INTEGER   AD1, AD2, ADDRSS, ADSAVE, COUNT, DAFILE, DUMMY, FREE, I, ID, ID1, IW(1), KFILE, LP, MAXL, MV, MXNEW, MXOLD,   &
                MXRST, N, NBADMU, NIL, NMXRST, NOACTN, NOR, NREAD, NWRITE, READID, READK, SAEVAL, SAVE, SAVFRE, SCPX,   &
                SOLCPX, V, WRITID, X, VER, ID2, LEN
LOGICAL   ALLOC, F, IO, MEQI, OPM, SAVE1, T, OPINV
COMMON   /STLMAC/ NOACTN, FREE, SAVE, SAVFRE
COMMON   /STLMAD/ ADDRSS, DAFILE, KFILE, LP, MAXL, NREAD, NWRITE
COMMON   /STLMID/ NIL, MV, V, MXNEW, MXOLD, MXRST, SCPX, SOLCPX, X
COMMON   /STLMIO/ SAEVAL, READID, WRITID, READK, N
COMMON   /STLMMI/ MEQI
COMMON   /STLMST/ TIME(24), COUNT(24), NBADMU, NMXRST, DUMMY
COMMON   /STLMTF/ T, F
COMMON   /STLMVR/ SMALL, VER
!
DATA      NOR  / 16 /
!
ST = SECOND()
COUNT(NOR) = COUNT(NOR) + 1
RANDVC = T
!
ID1 = ID
ID2 = V + 2
IF(.NOT. MEQI) GOTO 5
  ID2 = ID
  ID1 = V + 2
!
5     IF(.NOT. VER <= 2) GOTO 20
!     *************
!     VER = 1 OR 2.
!     *************
  IF(.NOT. ALLOC(ID1, NIL, AD1, AD2, NOACTN, W, IW(ADDRSS) ))GOTO 8888
!
  ADSAVE = AD1
  AD1    = AD1 - 1
!
  DO 10 I = 1, N
    AD1 = AD1 + 1
!         ******************************************
!         NOTE. RANFX SHOULD BE SUPPLIED BY THE USER.
!         ******************************************
    W(AD1) = RANFX(DUMMY)
10        CONTINUE
!
  GOTO 30
!
!     ********
!     VER = 3.
!     ********
20    CALL RAN3(ID1, N, W, IW, LEN)
IF(LEN /= 0) GOTO 8888
!
30    IF(.NOT. MEQI) GOTO 40
!       *****************
!       TRANSFORM VECTOR.
!       *****************
  IF(.NOT. OPINV(ID2, ID1, W, IW)) GOTO 8888
  CALL FREEID(ID1)
  GOTO 9999
!
!     *****************
!     TRANSFORM VECTOR.
!     *****************
40    IF(.NOT. OPM(ID2, ID1, W, IW)) GOTO 8888
IF(.NOT. OPINV(ID1, ID2, W, IW)) GOTO 8888
CALL FREEID(ID2)
!
IF(.NOT. VER <= 2) GOTO 9999
  IF(.NOT. SAVE1) GOTO 9999
!         ************
!         SAVE VECTOR.
!         ************
    IF(.NOT. IO(W(ADSAVE), ID, WRITID, N, IW(ADDRSS) )) GOTO 8888
!
    GOTO 9999
!
8888  IF(.NOT. VER == 3) GOTO 8889
  CALL ERROR(104, LEN)
  CALL ERROR(104, 104)
!
8889  RANDVC = F
CALL ERROR(NOR, NOR)
!
9999  TIME(NOR) = TIME(NOR) + SECOND() - ST
RETURN
END
