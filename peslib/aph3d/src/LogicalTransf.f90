LOGICAL FUNCTION LogicalTransf(ADD, SCOL, XORMX, VORMV, W, IW)
!
!   $RCSfile: LogicalTransf.f,v $ $Revision: 1.3 $
!   $Date: 89/07/28 09:57:48 $
!   $State: Stable $
!
!
! **********************************************************************
!
!     PURPOSE -
!
!         FOR A GIVEN S-VECTOR THIS ROUTINE COMPUTES X=V*S-ADD*V(P+1)
!         OR MX=MV*S-ADD*MV(P+1), WHERE ADD=-S(P)*BETA(P)/NU,
!         (ADD=-S(P)/NU IF ZERBET IS TRUE). V AND MV ARE MATRICES, SO
!         THE ROUTINE PERFORMS ESSENTIALLY A MATRIX VECTOR
!         MULTIPLICATION.
!
!
!     INPUT PARAMETERS -
!
!         ADD    = -S(P,I)*BETA(P)/NU(I) FOR AN I, IF ZERBET IS FALSE.
!                = -S(P, I)/NU(I), IF ZERBET = TRUE.
!         SCOL   = (S(1,I), ..., S(P,I))(T), I.E. COLUMN NUMBER I IN S.
!         XORMX  = THE COMPUTED X (OR MX) VECTOR GETS THE IDENTIFIER
!                  XORMX+1.
!         VORMV  = THE COLUMNS IN THE V (OR MV) MATRIX HAVE THE
!                  IDENTIFIERS, VORMV+1, ..., VORM+P, (VORMV+P+1).
!
!     PLEASE SEE THE PROGRAMMERS GUIDE FOR INFORMATION ABOUT
!     PARAMETERS NOT EXPLAINED ABOVE, AND FOR MORE DETAILS ABOUT
!     THE FUNCTION OF THE ROUTINE.
!
!
! **********************************************************************
!
USE Numeric_Kinds_Module
REAL(Kind=WP_Kind) ADD, SCOL(1), SECOND, ST, TIME, W(1)
INTEGER   CNEG, COUNT, CPOS, DUMMY, FREE, I, ITERNO,   IW(1), N, NBADMU, NMXRST, NOACTN, NOR, OLCPOS, P, REST,   &
                RNEW, ROLD, SAVE, SAVFRE, TCONV, VORMV, XORMX
LOGICAL   F, MULVEC, SUBVEC, T, USEMX, ZERBET
COMMON   /STLMAC/ NOACTN, FREE, SAVE, SAVFRE
COMMON   /STLMCT/ N, ITERNO, TCONV, CNEG, CPOS, OLCPOS, RNEW, ROLD, REST, P, USEMX, ZERBET
COMMON   /STLMST/ TIME(24), COUNT(24), NBADMU, NMXRST, DUMMY
COMMON   /STLMTF/ T, F
!
DATA      NOR  / 22 /
!
ST = SECOND()
COUNT(NOR) = COUNT(NOR) + 1
LogicalTransf = T
!     ************************
!     X = V1 * S1 (OR MX, MV).
!     ************************
IF(.NOT. MULVEC(XORMX, VORMV + 1, SCOL(1), W, IW, FREE, 1))GOTO 8888
!
IF(.NOT. P > 1) GOTO 50
  DO 40 I = 2, P
!         ********************************************
!         X = X + V2 * S2 + ... + VP * SP (OR MX, MV).
!         ********************************************
    IF(.NOT. SUBVEC(XORMX, VORMV + I, - SCOL(I), W, IW, FREE, 1))GOTO 8888
40        CONTINUE
!
!     *********************************
!     X = X - ADD * V(P+1) (OR MX, MV).
!     *********************************
50    IF(.NOT. SUBVEC(XORMX, VORMV + P + 1, ADD, W, IW, SAVFRE, 1))GOTO 8888
CALL FREEID(VORMV + P + 1)
GOTO 9999
!
8888  LogicalTransf = F
CALL ERROR(NOR, NOR)
!
9999  TIME(NOR) = TIME(NOR) + SECOND() - ST
RETURN
END
