LOGICAL FUNCTION SAVEXL(BETAPI, LAMBDA, NU, POINTR, S, MXNEG,  PMAX, CPALSO, W, IW)
USE Numeric_Kinds_Module
!
!   $RCSfile: savexl.f,v $ $Revision: 1.3 $
!   $Date: 89/07/28 09:57:37 $
!   $State: Stable $
!
!
! **********************************************************************
!
!     PURPOSE -
!
!         COMPUTES AND STORES EIGENVECTORS (AND MX-VECTORS IF USEMX).
!         SAVES EIGENVALUSES.
!
!
!     INPUT PARAMETERS -
!
!         MXNEG  = IF THE MX(CNEG) VECTORS ARE COMPUTED THEY GET
!                  IDENTIFIERS, MXNEG+1, ..., MXNEG+CNEG.
!         CPALSO = FALSE, COMPUTE AND SAVE CNEG-PAIRS. IF TRUE, HANDLE
!                  THE CPOS-PAIRS ALSO.
!
!     PLEASE SEE THE PROGRAMMERS GUIDE FOR INFORMATION ABOUT
!     PARAMETERS NOT EXPLAINED ABOVE, AND FOR MORE DETAILS ABOUT
!     THE FUNCTION OF THE ROUTINE.
!
!
! **********************************************************************
!
INTEGER   PMAX
REAL(Kind=WP_Kind) ADD, BETAPI(1), LAMBDA(1), NU(1), RDUMP, S(PMAX,1),  SECOND, ST, TIME, W(1)
INTEGER   ADDRSS, CNEG, COUNT, CPOS, DAFILE, DUMMY, ERRNO, I, I1,  I2, IDUMP, ITERNO, IW(1), J, J1, KFILE, LP, MAXL, MV,           &
                MXNEG, MXNEW, MXOLD, MXRST, MXS, N, NALSO, NBADMU, NIL,  NMXRST, NOR, NREAD, NWRITE, OLCPOS, P, POINT, POINTR(1),  &
                READID, READK, REST, RNEW, ROLD, SAEVAL, SCPX, SOLCPX, TCONV, V, WRITID, X
LOGICAL   CPALSO, F, T, LogicalTransf, USEMX, ZERBET
COMMON   /STLMAD/ ADDRSS, DAFILE, KFILE, LP, MAXL, NREAD, NWRITE
COMMON   /STLMCT/ N, ITERNO, TCONV, CNEG, CPOS, OLCPOS, RNEW, ROLD, REST, P, USEMX, ZERBET
COMMON   /STLMER/ RDUMP, ERRNO, IDUMP(2)
COMMON   /STLMID/ NIL, MV, V, MXNEW, MXOLD, MXRST, SCPX, SOLCPX, X
COMMON   /STLMIO/ SAEVAL, READID, WRITID, READK, NALSO
COMMON   /STLMST/ TIME(24), COUNT(24), NBADMU, NMXRST, DUMMY
COMMON   /STLMTF/ T, F
!
DATA      NOR  / 17 /
!
ST = SECOND()
COUNT(NOR) = COUNT(NOR) + 1
SAVEXL = T
!
!     ***********************
!     PREPARE FOR CNEG PAIRS.
!     ***********************
MXS = MXNEG
I1 = 1
I2 = CNEG
IF(CNEG == 0) GOTO 30
!
10    J = 1
!     *************************************
!     CHECK IF WE SHALL COMPUTE MX-VECTORS.
!     *************************************
IF(.NOT. USEMX) MXS = NIL
!
DO 20 I = I1, I2
  POINT = POINTR(I)
!       ********************
!       UPDATE TCONV.
!       COMPUTE EIGENVECTOR.
!       SAVE EIGENVALUE.
!       ********************
  TCONV = TCONV + 1
!       **********************************
!       CHECK IF ROOM TO STORE EIGENVALUE.
!       **********************************
  IF(.NOT. TCONV > MAXL) GOTO 15
    TCONV = TCONV - 1
    RDUMP = LAMBDA(POINT)
    CALL ERROR(NOR, 1)
    GOTO 8888
!
!       ****************
!       SAVE EIGENVALUE.
!       ****************
15      J1 = LP + TCONV - 1
  W(J1) = LAMBDA(POINT)
!
  ADD = -BETAPI(POINT)
  IF(ZERBET) ADD = -S(P, POINT) / NU(POINT)
  IF(.NOT. LogicalTransf(ADD, S(1, POINT), X + TCONV, V, W, IW)) GOTO 8888
!
  IF(.NOT. MXS /= NIL) GOTO 20
!         ******************
!         COMPUTE MX-VECTOR.
!         ******************
    IF(.NOT. LogicalTransf(ADD, S(1, POINT), MXS + J, MV, W, IW)) GOTO 8888
!
!
    J = J + 1
20      CONTINUE
!
!     ***************
!     CHECK IF READY.
!     ***************
30    IF(.NOT. (CPALSO .AND. I2 == CNEG .AND. CPOS > 0)) GOTO 9999
!       *************************
!       PREPARE FOR CPOS VECTORS.
!       *************************
  I1 = CNEG + 1
  I2 = CNEG + CPOS
  MXS = MXNEW
  GOTO 10
!
8888  SAVEXL = F
CALL ERROR(NOR, NOR)
!
9999  TIME(NOR) = TIME(NOR) + SECOND() - ST
RETURN
END
