LOGICAL FUNCTION IO(VEC, ID, CASE, NUMEL, ADDRSS)
USE Numeric_Kinds_Module
!
!   $RCSfile: io.f,v $ $Revision: 1.3 $
!   $Date: 89/07/28 09:57:17 $
!   $State: Stable $
!
!
! **********************************************************************
!
!     PURPOSE - (VER = 1 OR 2)
!
!         THIS ROUTINE READS OR WRITES A VECTOR OF LENGTH N
!         (VER = 1 OR 2).
!         IF VER = 1, IT ALSO READS THE K-MATRIX.
!
!     INPUT PARAMETERS -
!
!         VEC    = IS THE VECTOR TO BE WRITTEN.
!         ID     = IS THE IDENTIFIER OF THE VECTOR.
!         NUMEL  = THE LENGTH OF THE VECTOR.
!         CASE   = 1, NOT USED IN THE TOMS VERSION.
!                  2, READ THE VECTOR.
!                  3, WRITE THE VECTOR.
!                  4, READ THE K-MATRIX.
!         ADDRSS = ADDRESS VECTOR FOR THE N-VECTORS
!
!     OUTPUT PARAMETERS -
!
!         VEC    = IS THE VECTOR TO BE WRITTEN.
!
!     PLEASE SEE THE PROGRAMMERS GUIDE FOR INFORMATION ABOUT
!     PARAMETERS NOT EXPLAINED ABOVE, AND FOR MORE DETAILS ABOUT
!     THE FUNCTION OF THE ROUTINE.
!
!
! **********************************************************************
!
INTEGER   NUMEL
REAL(Kind=WP_Kind)  RDUMP, SECOND, ST, TIME, VEC(NUMEL)
INTEGER   ADDRS1, ADDRSS(1), CASE, COUNT, DAFILE, DUMMY, ERRNO,  ID, IDUMP, KFILE, LEFTP, LENADR, LP, MAXL, MAXREC, MV,      &
                MXNEW, MXOLD, MXRST, NBADMU, NIL, NMXRST, NOR, NREAD, NUMVEC, NWRITE, RIGHTC, RIGHTM, RIGHTP, SCPX, SOLCPX,  &
                STADEW, V, X
LOGICAL   F, T
COMMON   /STLMAD/ ADDRS1, DAFILE, KFILE, LP, MAXL, NREAD, NWRITE
COMMON   /STLMER/ RDUMP, ERRNO, IDUMP(2)
COMMON   /STLMEW/ LEFTP, LENADR, MAXREC, NUMVEC, RIGHTC, RIGHTM, RIGHTP, STADEW
COMMON   /STLMID/ NIL, MV, V, MXNEW, MXOLD, MXRST, SCPX, SOLCPX, X
COMMON   /STLMST/ TIME(24), COUNT(24), NBADMU, NMXRST, DUMMY
COMMON   /STLMTF/ T, F
DATA      NOR  / 9 /
!
ST = SECOND()
COUNT(NOR) = COUNT(NOR) + 1
IO = T
!
IF(CASE==2)THEN
!     *************
!     READ VEC(ID).
!     *************
100   IF(ADDRSS(ID) < 0) GOTO 9999
!     ********************************************************
!     THE IMPLEMENTOR SHOULD REPLACE THIS FORTRAN77 STATEMENT,
!     IF NECESSARY.
!     THE STATEMENT SHOULD READ VEC OF LENGTH NUMEL (=N), WITH
!     INDEX ADDRSS(ID) ON THE DIRECT ACCESS FILE DAFILE.
!     IN CASE OF ERROR, GOTO 110.
!     ********************************************************
READ(DAFILE, REC = ADDRSS(ID), ERR = 110) VEC
NREAD = NREAD + 1
GOTO 9999
!
110   IDUMP(1) = ADDRSS(ID)
CALL ERROR(NOR, 2)
GOTO 8888
!     **************
!     WRITE VEC(ID).
!     **************
ELSEIF(CASE==3)THEN
200   IF(ADDRSS(ID) < 0) GOTO 9999
IF(.NOT. ADDRSS(ID) == 0) GOTO 202
  IF(RIGHTP == RIGHTC) RIGHTP = RIGHTM
  IF(RIGHTP < LEFTP .OR. RIGHTM < LEFTP) GOTO 220
  IF(.NOT. ID <= X) GOTO 201
    ADDRSS(ID) = RIGHTP
    RIGHTP = RIGHTP - 1
    GOTO 202
!
201     ADDRSS(ID) = LEFTP
  LEFTP = LEFTP + 1
!
!     *****************************************************************
!     THE IMPLEMENTOR SHOULD REPLACE THE FOLLOWING FORTRAN77 STATEMENT,
!     IF NECESSARY.
!     THE STATEMENT SHOULD WRITE VEC OF LENGTH NUMEL (=N), WITH INDEX
!     ADDRSS(ID) ON THE DIRECT ACCESS FILE DAFILE.
!     IN CASE OF ERROR GOTO 210.
!     *****************************************************************
202   CONTINUE
WRITE(DAFILE, REC = ADDRSS(ID), ERR = 210) VEC
NWRITE = NWRITE + 1
GOTO 9999
!
210   IDUMP(1) = ADDRSS(ID)
CALL ERROR(NOR, 3)
GOTO 8888
220   CALL ERROR(NOR, 1)
GOTO 8888
!     **************
!     READ K-MATRIX.
!     **************
ELSEIF(CASE==4)THEN
300   REWIND KFILE
READ(KFILE) VEC
GOTO 9999
ENDIF
!
8888  IO = F
CALL ERROR(NOR, NOR)
!
9999  TIME(NOR) = TIME(NOR) + SECOND() - ST
RETURN
END
