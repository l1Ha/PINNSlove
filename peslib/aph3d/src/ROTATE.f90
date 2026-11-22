SUBROUTINE ROTATE(S, POINTR, FINAL, PMAX, NUMDEL)
USE Numeric_Kinds_Module
!
!   $RCSfile: rotate.f,v $ $Revision: 1.3 $
!   $Date: 89/07/28 09:57:37 $
!   $State: Stable $
!
!
! **********************************************************************
!
!     PURPOSE -
!
!         GIVEN TWO S-VECTORS CORRESPONDING TO NUMERICALLY INDENTICAL
!         NU-EIGENVALUES THIS ROUTINE CONSTRUCTS ONE GOOD VECTOR. THIS
!         IS MADE WITH A ROTATION MATRIX.
!
!
!     INPUT PARAMETERS -
!
!         POINTR = A PART OF THE POINTR VECTOR.
!         FINAL  = TRUE IF THE LAST LANCZOS STEP HAS BEEN TAKEN.
!
!
!     OUTPUT PARAMETERS -
!
!         NUMDEL = UPDATED VALUE ON THE NUMBER OF DELETED AND
!                  CONVERGED VECTORS.
!
!     PLEASE SEE THE PROGRAMMERS GUIDE FOR INFORMATION ABOUT
!     PARAMETERS NOT EXPLAINED ABOVE, AND FOR MORE DETAILS ABOUT
!     THE FUNCTION OF THE ROUTINE.
!
!
!
! **********************************************************************
!
INTEGER   PMAX
REAL(Kind=WP_Kind)  H1, H2, Q, S(PMAX, 2)
INTEGER   BAD, CNEG, CPOS, GOOD, I, ITERNO, N, NUMDEL,OLCPOS, P, POINTR(2), REST, RNEW, ROLD, STEP, TCONV
LOGICAL   FINAL, USEMX, ZERBET
COMMON   /STLMCT/ N, ITERNO, TCONV, CNEG, CPOS, OLCPOS, RNEW, ROLD, REST, P, USEMX, ZERBET
!
!     *****************************************************
!     IF( ... )THEN WE HAVE DELETED A CONVERGED EIGENPAIR.
!     *****************************************************
IF(POINTR(1) > 0 .AND. POINTR(2) > 0) NUMDEL = NUMDEL + 1
I = 1
!     *******************************************************
!     THE BETTER VECTOR (GOOD) HAS A GREATER FIRST COMPONENT.
!     *******************************************************
GOOD = IABS(POINTR(1))
BAD = IABS(POINTR(2))
IF(.NOT. ABS(S(1, GOOD)) < ABS(S(1, BAD))) GOTO 10
  I = GOOD
  GOOD = BAD
  BAD = I
  I = 2
10    POINTR(I) = GOOD
I = 3 - I
!     *********************
!     CLEAR THE BAD POINTR.
!     *********************
POINTR(I) = 0
!
!     ******************************
!     CONSTRUCT THE ROTATION MATRIX.
!     ******************************
Q = S(1, BAD) / S(1, GOOD)
H2 = SQRT(1.0D0  + Q * Q)
H1 = 1.0D0  / H2
H2 = Q / H2
!
!     ******************************************************
!     IF NOT FINAL, CHANGE ONLY THE TOP AND BOTTOM ELEMENTS.
!     ******************************************************
STEP = P - 1
IF(FINAL) STEP = 1
!
!     *******
!     ROTATE.
!     *******
DO 20 I = 1, P, STEP
  S(I, GOOD) = H1 * S(I, GOOD) + H2 * S(I, BAD)
20      CONTINUE
!
RETURN
END
