SUBROUTINE OPMSUB(MX, X, M, D, N, DIAGM)
USE Numeric_Kinds_Module
!
!   $RCSfile: opmsub.f,v $ $Revision: 1.3 $
!   $Date: 89/07/28 09:57:29 $
!   $State: Stable $
!
!
! **********************************************************************
!
!     PURPOSE - (VER = 1)
!
!         MX=M*X IS COMPUTED.
!
!
!     INPUT PARAMETERS -
!
!         D     = POINTER VECTOR TO DIAGONAL ELEMENTS OF M (IF
!                 M IS NOT DIAGONAL).
!         DIAGM = PROFIL <= 2
!         M     = M MATRIX STORED ACCORDING TO THE USER GUIDE.
!         N     = DIMENSION OF M MATRIX.
!         X     = VECTOR TO BE MULTIPLIED BY M.
!
!
!     OUTPUT PARAMETERS -
!
!         MX = RESULTING VECTOR.
!
! **********************************************************************
!
REAL(Kind=WP_Kind)  M(1), MX(1), SCPR, X(1)
INTEGER   D(1), DI, DIM1, I, LENM1, N, ROWNO, TOP
LOGICAL   DIAGM
!
IF(.NOT. DIAGM) GOTO 40
!       ***********************************
!       M IS DIAGONAL (BUT NOT EQUAL TO I).
!       ***********************************
  DO 30 I = 1, N
    MX(I) = M(I) * X(I)
30        CONTINUE
  GOTO 9999
!
!     ******************
!     M IS NON DIAGONAL.
!     ******************
40    DIM1 = 0
!
DO 50 I = 1, N
  TOP = DIM1 + 1
  DI = D(I)
  LENM1 = DI - TOP
  ROWNO = I - LENM1
  MX(I) = SCPR(X(ROWNO), M(TOP), LENM1 + 1)
  IF(LENM1 > 0) CALL SUBV(MX(ROWNO), M(TOP), - X(I), LENM1)
  DIM1 = DI
50      CONTINUE
!
9999  RETURN
END
