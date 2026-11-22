SUBROUTINE SPL1D2 (N, X, F, W, IJ, Y, TAB)
USE Numeric_Kinds_Module
!***BEGIN PROLOGUE  SPL1D2
!***DATE WRITTEN   790301   (YYMMDD)
!***REVISION DATE  831207   (YYMMDD)
!***CATEGORY NO.  E1A
!***KEYWORDS  SPLINES,CUBIC SPLINES,INTERPOLATION,DERIVATIVES
!***AUTHOR  JORDAN, TOM, LOS ALAMOS NATIONAL LABORATORY
!***PURPOSE  Interpolates values of a function and its first and
!            second derivatives.
!***DESCRIPTION
!
!     SPL1D2 interpolates in a table, returning values of the function
!     and its first and second derivatives.  This is normally used
!     in conjunction with SPL1D1 which will supply second derivatives.
!
!     -ARGUMENTS-
!
!      N = INTEGER specifying the number of points.
!      X = origin of the table of floating point values of independent
!          variable.  This table must be in ascending order.
!      F = origin of the table of floating point values of the dependent
!          variable.
!      W = origin of the table of floating point values of the second
!          derivative.  (Normally supplied through use of SPL1D1.)
!     IJ = INTEGER specifying spacing in tables F and W.
!      Y = floating point value at which interpolation is desired.
!    TAB = array of dimension >= 3 which contains the results of the
!          interpolation.
!
!     Upon return from SPL1D2, TAB(1) contains F(Y), TAB(2) contains
!     F'(Y), and TAB(3) contains F''(Y).  There are no error messages
!     in SPL1D2.
!
!     NOTE: If Y is outside the range of the X table, extrapolation
!     occurs.
!
!     METHOD:
!
!     A complete discussion of the method may be found in the reference
!     listed below.
!
!***REFERENCES  Walsh, Ahlberg, and Nilson, Journal of Mathematics and
!                 Mechanics, Vol. II, No. 2, 1962.
!***ROUTINES CALLED  (NONE)
!***END PROLOGUE  SPL1D2
IMPLICIT NONE
INTEGER I, IJ, K1, MI, N, ilo
REAL(Kind=WP_Kind) A, B, C, FLK, Y
REAL(Kind=WP_Kind) X (N), F (N), W (N), TAB (3)
SAVE ilo
!***FIRST EXECUTABLE STATEMENT  SPL1D2
!
!     LOCATE Y IN THE X TABLE
!
DATA ilo / 1 /
!
CALL FndInt (N - 1, X, Y, ilo, I)
I = MAX0 (I, 1)
MI = (I - 1) * IJ + 1
K1 = MI + IJ
FLK = X (I + 1) - X (I)
!
!     CALCULATE F(Y)
!
A = (W (MI) * (X (I + 1) - Y) **3 + W (K1) * (Y - X (I) ) **3) &
 / (6. * FLK)
B = (F (K1) / FLK - W (K1) * FLK / 6.) * (Y - X (I) )
C = (F (MI) / FLK - FLK * W (MI) / 6.) * (X (I + 1) - Y)
TAB (1) = A + B + C
!
!     CALCULATE THE FIRST DERIVATIVE AT Y
!
A = (W (K1) * (Y - X (I) ) **2 - W (MI) * (X (I + 1) - Y) **2) &
 / (2. * FLK)
B = (F (K1) - F (MI) ) / FLK
C = FLK * (W (MI) - W (K1) ) / 6.
TAB (2) = A + B + C
!
!     CALCULATE THE SECOND DERIVATIVE AT Y
!
TAB (3) = (W (MI) * (X (I + 1) - Y) + W (K1) * (Y - X (I) ) ) &
 / FLK
RETURN
ENDSUBROUTINE SPL1D2
