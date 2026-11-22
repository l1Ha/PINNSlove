LOGICAL FUNCTION SUBVEC(ID1, ID2, C, W, IW, ACTION, CASE)
USE Numeric_Kinds_Module
!
!   $RCSfile: subvec.f,v $ $Revision: 1.3 $
!   $Date: 89/07/28 09:57:47 $
!   $State: Stable $
!
!
! **********************************************************************
!
!     PURPOSE -
!
!         DRIVER ROUTINE FOR SUBV. COMPUTES
!         VEC(ID1)=VEC(ID1)-C*VEC(ID2).
!
!
!     INPUT PARAMETERS -
!
!         ID1    = IDENTIFIER OF THE RESULTING VECTOR.
!         ID2    = IDENTIFIER OF THE INCOMING VECTOR.
!         C      = REAL CONSTANT.
!         ACTION = DETERMINES IF VEC(ID1) OR VEC(ID2) SHOULD
!                  BE SAVED OR NOT.
!         CASE   = NOT USED IN THIS VERSION.
!
!     PLEASE SEE THE PROGRAMMERS GUIDE FOR INFORMATION ABOUT
!     PARAMETERS NOT EXPLAINED ABOVE, AND FOR MORE DETAILS ABOUT
!     THE FUNCTION OF THE ROUTINE.
!
!
! **********************************************************************
!
REAL(Kind=WP_Kind)   C, SECOND, ST, TIME, W(1), SMALL
INTEGER   ACTION, AD1, AD2, ADDRSS, CASE, COUNT, DAFILE, DUMMY, FREE, ID1, ID2, IW(1), KFILE, LP, MAXL, N, NBADMU,   &
                NMXRST, NOACTN, NOR, NREAD, NWRITE, READID, READK, SAEVAL, SAVE, SAVFRE, WRITID, VER, LEN
LOGICAL   ALLOC, F, IO, T
COMMON   /STLMAC/ NOACTN, FREE, SAVE, SAVFRE
COMMON   /STLMAD/ ADDRSS, DAFILE, KFILE, LP, MAXL, NREAD, NWRITE
COMMON   /STLMIO/ SAEVAL, READID, WRITID, READK, N
COMMON   /STLMST/ TIME(24), COUNT(24), NBADMU, NMXRST, DUMMY
COMMON   /STLMTF/ T, F
COMMON   /STLMVR/ SMALL, VER
!
DATA      NOR  / 21 /
!
ST = SECOND()
COUNT(NOR) = COUNT(NOR) + 1
SUBVEC = T
!
IF(.NOT. VER <= 2) GOTO 10
!       *************
!       VER = 1 OR 2.
!       *************
  IF(.NOT. ALLOC(ID1, ID2, AD1, AD2, ACTION, W, IW(ADDRSS) )) GOTO 8888
!
  CALL SUBV(W(AD1), W(AD2), C, N)
!
  IF(.NOT. (ACTION == SAVE .OR. ACTION == SAVFRE)) GOTO 9999
    IF(.NOT. IO(W(AD1), ID1, WRITID, N, IW(ADDRSS) )) GOTO 8888
    GOTO 9999
!
!     ********
!     VER = 3.
!     ********
10    CALL SUB3(ID1, ID2, C, N, W, IW, LEN)
!
IF(.NOT. LEN /= 0) GOTO 9999
  CALL ERROR(106, LEN)
  CALL ERROR(106, 106)
!
8888  SUBVEC = F
CALL ERROR(NOR, NOR)
!
9999  TIME(NOR) = TIME(NOR) + SECOND() - ST
RETURN
END
