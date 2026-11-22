LOGICAL FUNCTION PREPSV(IDX, IDMX, NUMBER, W, IW)
USE Numeric_Kinds_Module
!
!   $RCSfile: prepsv.f,v $ $Revision: 1.3 $
!   $Date: 89/07/28 09:57:30 $
!   $State: Stable $
!
!
! **********************************************************************
!
!     PURPOSE -
!
!         PREPARES THE STARTINGVECTOR, I.E. ORTHOGONALIZES IT AGAINST
!         SOME ALREADY COMPUTED EIGENVECTORS.
!
!
!
!     INPUT PARAMETERS -
!
!         IDX    = IDX+1, ..., IDX+NUMBER ARE THE IDENTIFIERS FOR THE
!                  X-VECTORS.
!         IDMX   = IDMX+1, ..., IDMX+NUMBER ARE THE IDENTIFIERS
!                  FOR THE MX-VECTORS (IF USED).
!         NUMBER = IS THE NUMBER OF VECTORS.
!
!     PLEASE SEE THE PROGRAMMERS GUIDE FOR INFORMATION ABOUT
!     PARAMETERS NOT EXPLAINED ABOVE, AND FOR MORE DETAILS ABOUT
!     THE FUNCTION OF THE ROUTINE.
!
!
! **********************************************************************
!
REAL(Kind=WP_Kind)  DOT, SECOND, ST, TIME, W(1)
INTEGER   CNEG, COUNT, CPOS, DUMMY, FREE, I, IDMX, IDX, ITERNO,     &
                IW(1), MV, MXNEW, MXOLD, MXRST, N, NBADMU, NIL, NMXRST, NOACTN, NOR, NUMBER, OLCPOS, P, REST, RNEW, ROLD, SAVE,   &
                SAVFRE, SCPX, SOLCPX, TCONV, V, X
LOGICAL   F, MEQI, OPM, SCPROD, SUBVEC, T, USEMX, ZERBET
COMMON   /STLMAC/ NOACTN, FREE, SAVE, SAVFRE
COMMON   /STLMCT/ N, ITERNO, TCONV, CNEG, CPOS, OLCPOS, RNEW, ROLD, REST, P, USEMX, ZERBET
COMMON   /STLMID/ NIL, MV, V, MXNEW, MXOLD, MXRST, SCPX, SOLCPX, X
COMMON   /STLMMI/ MEQI
COMMON   /STLMST/ TIME(24), COUNT(24), NBADMU, NMXRST, DUMMY
COMMON   /STLMTF/ T, F
!
DATA      NOR  / 15 /
!
ST = SECOND()
COUNT(NOR) = COUNT(NOR) + 1
PREPSV = T
!
IF(NUMBER <= 0) GOTO 9999
!
IF(.NOT. USEMX) GOTO 15
!       ************************
!       THE MX-VECTORS ARE USED.
!       ************************
  DO 10 I = 1, NUMBER
!         *****************************
!         V1 = V1 - (V1(T) * (MX)) * X.
!         *****************************
    IF(.NOT. SCPROD(V + 1, IDMX + I, DOT, W, IW, FREE, 2))GOTO 8888
!
    IF(.NOT. SUBVEC(V + 1, IDX + I, DOT, W, IW, FREE, 3))GOTO 8888
!
10        CONTINUE
  GOTO 9999
!
15    IF(.NOT. MEQI) GOTO 30
!       ******
!       M = I.
!       ******
  DO 20 I = 1, NUMBER
!         **************************
!         V1 = V1 - (V1(T) * X) * X.
!         **************************
    IF(.NOT. SCPROD(V + 1, IDX + I, DOT, W, IW, NOACTN, 2))GOTO 8888
!
    IF(.NOT. SUBVEC(V + 1, IDX + I, DOT, W, IW, FREE, 3))GOTO 8888
!
20      CONTINUE
GOTO 9999
!
!     ********************
!     M IS NOT EQUAL TO I.
!     ********************
30    CALL FREEID(V + 1)
DO 40 I = 1, NUMBER
!       ****************************************************
!       V2 = M * X, (V2 = V + 2 IS USED AS A SCRATCH VECTOR)
!       V1 = V1 - (V1(T) * V2) * X.
!       ****************************************************
  IF(.NOT. OPM(V + 2, IDX + I, W, IW)) GOTO 8888
  CALL FREEID(IDX + I)
  IF(.NOT. SCPROD(V + 2, V + 1, DOT, W, IW, NOACTN, 2)) GOTO 8888
  CALL FREEID(V + 2)
  IF(.NOT. SUBVEC(V + 1, IDX + I, DOT, W, IW, SAVFRE, 2))GOTO 8888
!
  CALL FREEID(IDX + I)
40      CONTINUE
GOTO 9999
!
8888  PREPSV = F
CALL ERROR(NOR, NOR)
!
9999  TIME(NOR) = TIME(NOR) + SECOND() - ST
RETURN
END
