LOGICAL FUNCTION ALLOC(ID1, ID2, AD1, AD2, ACTION, W, ADDRSS)
USE Numeric_Kinds_Module
SAVE
!
!   $RCSfile: alloc.f,v $ $Revision: 1.3 $
!   $Date: 89/07/28 09:57:03 $
!   $State: Stable $
!
!
! **********************************************************************
!
!     PURPOSE - (VER = 1 OR 2)
!
!         THIS ROUTINE ESTABLISHES CONNECTIONS BETWEEN IDENTIFIERS
!         AND VECTORS OF LENGTH N IN W.
!         IF VEC(ID2) DOES NOT LIE IN W ALREADY IT IS
!         FETCHED FROM SECONDARY STORAGE WITH THE HELP OF IO. IF AD IS
!         THE ADDRESS IN W THEN
!         VEC(ID)=(W(AD), ..., W(AD+N-1)), FOR ID = ID1, ID2, AND
!         AD = AD1, AD2.
!
!     INPUT PARAMETERS -
!
!         ID1, ID2   IDENTIFIERS OF THE VECTORS
!         ACTION     THE CORRESPONDING ACTION
!         ADDRSS     ADDRESS VECTOR
!
!
!     OUTPUT PARAMETERS -
!
!         AD1, AD2   THE ADDRESSES CORRESPONDING TO ID1 AND ID2.
!
!
!     PLEASE SEE THE PROGRAMMERS GUIDE FOR INFORMATION ABOUT
!     PARAMETERS NOT EXPLAINED ABOVE, AND FOR MORE DETAILS ABOUT
!     THE FUNCTION OF THE ROUTINE.
!
!
! **********************************************************************
!
REAL(Kind=WP_Kind) RDUMP, SECOND, ST, TIME, W(1)
INTEGER   ACTION, ACTIVE, AD(2), AD1, AD2, ADDRSS(1), COUNT, DUMMY, ERRNO, FREE, FREPOS, HIT, I, ID, ID1, ID2, IDUMP,       &
          MV, MXNEW, MXOLD, MXRST, N, NBADMU, NIL, NMXRST, NOACTN,  NOR, READID, READK, SAEVAL, SAVE, SAVFRE, SCPX, SOLCPX, &
          V, WAD, WRITID, X
LOGICAL   F, IO, T
COMMON   /STLMAC/ NOACTN, FREE, SAVE, SAVFRE
COMMON   /STLMER/ RDUMP, ERRNO, IDUMP(2)
COMMON   /STLMID/ NIL, MV, V, MXNEW, MXOLD, MXRST, SCPX, SOLCPX, X
COMMON   /STLMIO/ SAEVAL, READID, WRITID, READK, N
COMMON   /STLMST/ TIME(24), COUNT(24), NBADMU, NMXRST, DUMMY
COMMON   /STLMTF/ T, F
COMMON   /STLMWH/ WAD(2), ACTIVE(2)
!
DATA      NOR  / 1 /
!
ST = SECOND()
COUNT(NOR) = COUNT(NOR) + 1
ALLOC = T
!
AD(1) = NIL
AD(2) = NIL
ID = ID1
!     *********************
!     LOOP FOR ID1 AND ID2.
!     *********************
I = 1
!
10    FREPOS = 0
IF(ACTIVE(1) == NIL) FREPOS = 1
IF(ACTIVE(2) == NIL) FREPOS = 2
!
!     ******************************
!     CHECK IF ID ALREADY LIES IN W.
!     ******************************
HIT = 0
IF(ID == ACTIVE(1)) HIT = 1
IF(ID == ACTIVE(2)) HIT = 2
!
IF(.NOT. ADDRSS(ID) < 0) GOTO 15
  AD(I) = -ADDRSS(ID)
  GOTO 40
!
15    IF(.NOT. HIT /= 0) GOTO 20
  AD(I) = WAD(HIT)
  GOTO 40
!
20      IF(.NOT. FREPOS <= 0) GOTO 30
!         ******************************************
!         NO ROOM LEFT. SHOULD NEVER HAPPEN, BUT ...
!         ******************************************
!
    IDUMP(1) = ID
    IDUMP(2) = ACTION
    CALL ERROR(NOR, 1)
    GOTO 8888
!
!     *******
!     UPDATE.
!     *******
30    ACTIVE(FREPOS) = ID
AD(I) = WAD(FREPOS)
IF(.NOT. I == 2) GOTO 40
!       ************
!       GET VEC(ID).
!       ************
  AD2 = AD(2)
  IF(.NOT. IO(W(AD2), ID, READID, N, ADDRSS)) GOTO 8888
!
40    IF(I == 2) GOTO 50
I = 2
ID = ID2
IF(ID == NIL) GOTO 50
GOTO 10
!
50    IF(ACTION == NOACTN .OR. ACTION == SAVE) GOTO 9999
!
!       **********
!       DELETE ID.
!       **********
  IF(ACTION == SAVFRE) ID = ID1
  IF(ID == ACTIVE(1)) ACTIVE(1) = NIL
  IF(ID == ACTIVE(2)) ACTIVE(2) = NIL
!
GOTO 9999
!
8888  ALLOC = F
CALL ERROR(NOR, NOR)
!
9999  TIME(NOR) = TIME(NOR) + SECOND() - ST
AD1 = AD(1)
AD2 = AD(2)
!
RETURN
END FUNCTION ALLOC
