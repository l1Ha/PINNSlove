SUBROUTINE NEWMU(LAMBDA, POINTR, MAXNU, LGEB)
USE Numeric_Kinds_Module
!
!   $RCSfile: newmu.f,v $ $Revision: 1.3 $
!   $Date: 89/07/28 09:57:25 $
!   $State: Stable $
!
!
! **********************************************************************
!
!     PURPOSE -
!
!         COMPUTES THE NEXT SHIFT. THERE ARE THREE CASES.
!
!         EIGENVALUES >= B HAVE CONVERGED. SET NEXT SHIFT TO B.
!
!         WE HAVE AT LEAST ONE CONVERGED EIGENVALUE GREATER
!         THE SHIFT.
!
!         WE HAVE ONLY NOT CONVERGED EIGENVALUES GREATER THAN THE SHIFT.
!
!
!     INPUT PARAMETERS -
!
!         LGEB  = T IF A CONVERGED EIGENVALUE >= B HAS BEEN COMPUTED AND
!                 THROWN AWAY. FALSE OTHERWISE.
!         MAXNU = MAX(NU(I)), I=1, ..., P.
!
!     PLEASE SEE THE PROGRAMMERS GUIDE FOR INFORMATION ABOUT
!     PARAMETERS NOT EXPLAINED ABOVE, AND FOR MORE DETAILS ABOUT
!     THE FUNCTION OF THE ROUTINE.
!
!
! **********************************************************************
!
REAL(Kind=WP_Kind)   ALTMU, LAMBDA(1), MAXNU, MU, NEXTMU, OLDMU, A, B
INTEGER   CNEG, CPOS, ITERNO, N, OLCPOS, P, POINT, POINTR(1),  REST, RNEW, ROLD, TCONV, NUMEIG, MAXW, MAXIW
LOGICAL   USEMX, ZERBET, T, F, REACHB, USSMXR, USEDB, LGEB
COMMON   /STLMCT/ N, ITERNO, TCONV, CNEG, CPOS, OLCPOS, RNEW, ROLD, REST, P, USEMX, ZERBET
COMMON   /STLMIN/ A, B, NUMEIG, MAXW, MAXIW
COMMON   /STLMMU/ MU, OLDMU, NEXTMU, ALTMU
COMMON   /STLMTF/ T, F
COMMON   /STLMTS/ REACHB, USSMXR, USEDB
!
!     ***************************************************
!     CONVERGED EIGENVALUES GREATER THAN B HAVE OCCURRED.
!     SET THE NEXT SHIFT TO B.
!     ***************************************************
IF(.NOT. LGEB) GOTO 5
  REACHB = T
  NEXTMU = B
  ALTMU  = MU + 0.1D0 * (B - MU)
  GOTO 9999
!
5     IF(.NOT. CPOS > 0) GOTO 10
!       ***********
!       FIRST CASE.
!       ***********
  POINT = CNEG + CPOS
  POINT = POINTR(POINT)
!       *****************************************************
!       THE GREATEST LAMBDA IS BISECTING THE DISTANCE BETWEEN
!       THE NEXT SHIFT (NEXTMU) AND MU.
!       ALTMU IS USED IN CONNECTION WITH MXREST.
!       *****************************************************
!       NEXTMU = 2.0D0 * LAMBDA(POINT) - MU                              .
!       ALTMU = 1.1D0 * LAMBDA(POINT) - 0.1D0 * MU                       .
  NEXTMU = 1.2D0 * LAMBDA(POINT) - 0.2D0 * MU
  ALTMU = 1.05D0 * LAMBDA(POINT) - 0.05D0 * MU
!       ***************************
!       CHECK IF WE HAVE REACHED B.
!       ***************************
  IF(.NOT. NEXTMU >= B) GOTO 9999
    REACHB = T
    NEXTMU = B
    GOTO 9999
!
!     ************
!     SECOND CASE.
!     ************
10    NEXTMU = MU + 0.9D0 / MAXNU
!     ***************************
!     CHECK IF WE HAVE REACHED B.
!     ***************************
IF(.NOT. NEXTMU >= B) GOTO 20
  REACHB = T
  NEXTMU = B
20    ALTMU = 0.9D0 * MU + 0.1D0 * NEXTMU
!
9999  RETURN
END
