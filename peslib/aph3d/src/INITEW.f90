SUBROUTINE INITEW(ADDRSS)
USE Numeric_Kinds_Module
!
!   $RCSfile: initew.f,v $ $Revision: 1.3 $
!   $Date: 89/07/28 09:57:15 $
!   $State: Stable $
!
!
! **********************************************************************
!
!     PURPOSE - (VER = 1 OR 2)
!
!         THIS ROUTINE UPDATES THE ADDRESS VECTOR TO THE
!         DIRECT ACCESS FILE.
!
!     PLEASE SEE THE PROGRAMMERS GUIDE FOR INFORMATION ABOUT
!     PARAMETERS NOT EXPLAINED ABOVE, AND FOR MORE DETAILS ABOUT
!     THE FUNCTION OF THE ROUTINE.
!
!
! **********************************************************************
!
REAL(Kind=WP_Kind)  FC, FN
INTEGER   ADDRSS(1), ADR, CNEG, COPT, CPOS, I, I1, ITERNO, LEFTP,  LENADR, MAXREC, MV, MXNEW, MXOLD, MXREST, MXRST, N, NIL,   &
                NUMMV, NUMMXR, NUMV, NUMVEC, OLCPOS, P, PMAX, POPT, REMAIN, REST, RIGHTC, RIGHTM, RIGHTP, RNEW, ROLD, SCPX,   &
                SOLCPX, STADEW, TCONV, V, X
LOGICAL   MEQI, USEMX, ZERBET
COMMON   /STLMCT/ N, ITERNO, TCONV, CNEG, CPOS, OLCPOS, RNEW, ROLD, REST, P, USEMX, ZERBET
COMMON   /STLMEW/ LEFTP, LENADR, MAXREC, NUMVEC, RIGHTC, RIGHTM, RIGHTP, STADEW
COMMON   /STLMID/ NIL, MV, V, MXNEW, MXOLD, MXRST, SCPX, SOLCPX, X
COMMON   /STLMMI/ MEQI
COMMON   /STLMPL/ PMAX, POPT, COPT, MXREST
!
!
IF(.NOT. USEMX) GOTO 25
!       *****************************
!       SET IF OLCPOS SHOULD BE ZERO.
!       *****************************
  RIGHTC = MAXREC + 1
  RIGHTM = RIGHTC
  IF(.NOT. OLCPOS > 0) GOTO 25
!
!         **************************************
!         ERASE ALL BUT THE X AND MXOLD VECTORS.
!         **************************************
    DO 21 I = 1, MXOLD
      ADDRSS(I) = 0
21          CONTINUE
!
    I1 = MXOLD + OLCPOS + 1
    DO 22 I = I1, X
      ADDRSS(I) = 0
22          CONTINUE
!
!         ***********************************************
!         RIGHTC POINTS TO THE  RECORD HOLDING  THE FIRST
!         VECTOR IN MXOLD. RIGHTM POINTS TO THE LEFT OF
!         THE LAST VECTOR IN MXOLD.
!         ***********************************************
    RIGHTM = ADDRSS(I1 - 1) - 1
    RIGHTC = ADDRSS(MXOLD + 1)
!
    GOTO 35
!
!
!     ************************************
!     WE COME HERE IF USEMX = F. ERASE ALL
!     BUT THE X-VECTORS.
!     ************************************
25    DO 30 I = 1, X
  ADDRSS(I) = 0
30      CONTINUE
!     *************
!     RESET RIGHTP.
!     *************
35    RIGHTP = MAXREC
!
!     *************************************
!     CHECK IF EVERYTHING LIES ON THE FILE.
!     *************************************
IF(NUMVEC == 0) GOTO 9999
!
!     *********************************************
!     MAKE AN OPTIMAL DIVISION OF THE STORAGE IN W.
!     *********************************************
FN = NUMVEC
FC = COPT
IF(MEQI) NUMMV = 0
IF(.NOT. MEQI .AND. .NOT. USEMX) NUMMV = (FN - FC) * 0.5D0
IF(.NOT. MEQI .AND.       USEMX) NUMMV = (FN - FC * 0.5D0) * 0.5D0
!
IF(NUMMV < 0) NUMMV = 0
!     ************************************
!     COMPUTE OPTIMAL NUMBER OF V-VECTORS.
!     CHECK IF GREATER THAN PMAX+1.
!     ************************************
NUMV = NUMVEC - NUMMV
I = PMAX + 1
IF(NUMV > I) NUMV = I
IF(NUMMV > I) NUMMV = I
!
REMAIN = NUMVEC - (NUMV + NUMMV)
IF(.NOT. USEMX) GOTO 50
!       ***********************************************
!       IF THERE REMAIN VECTORS WE PLACE THEM IN MXRST,
!       IF USEMX = T.
!       ***********************************************
  IF(REMAIN > PMAX) REMAIN = PMAX
  NUMMXR = REMAIN
!
!
50    ADR = - STADEW
!     ***************************
!     HERE ADDRSS IS INITIALIZED.
!     ***************************
CALL INITAD(NUMMV, MV, ADR, N, ADDRSS)
CALL INITAD(NUMV, V, ADR, N, ADDRSS)
IF(USEMX) CALL INITAD(NUMMXR, MXRST, ADR, N, ADDRSS)
!
!
9999  RETURN
END
