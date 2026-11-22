SUBROUTINE SPL1D1 (N, X, F, W, IOP, IJ, A, B, C)
USE Numeric_Kinds_Module
SAVE
!***BEGIN PROLOGUE  SPL1D1
!***DATE WRITTEN   790301   (YYMMDD)
!***REVISION DATE  831207   (YYMMDD)
!***CATEGORY NO.  E1A
!***KEYWORDS  SPLINES,CUBIC SPLINES,COEFFICIENTS
!***AUTHOR  JORDAN, TOM, LOS ALAMOS NATIONAL LABORATORY
!***PURPOSE  Finds coefficients of cubic for spline interpolation.
!***DESCRIPTION
!
!     SPL1D1 computes the coefficients of a sequence of cubics for
!     spline interpolation.  This SUBROUTINE will normally be used in
!     conjunction with SPL1D2.
!
!     -ARGUMENTS-
!
!      N = INTEGER specifying the number of points.
!      X = origin of the table of floating point values of independent
!          variable.  This table must be in ascending order.
!      F = origin of the table of floating point values of the dependent
!          variable.
!      W = origin of the table of floating point values of the second
!          derivatives of the spline at each knot.  The W's and F's
!          define the cubic over each interval.
!    IOP = an array of dimension 2 containing combinations of the
!          integers 1 through 5 for specifying the boundary conditions
!          (see Boundary Conditions).
!     IJ = INTEGER specifying spacing in tables F and W.  Normally,
!          IJ = 1.
!  A,B,C = arrays of dimension N used for temporary storage.  These
!          arrays may be used by the calling program at any other
!          time.
!
!     The second derivatives are computed and stored in W, length N.
!
!     ERROR CONDITIONS:
!
!        If option 5 is specified with N < 4, this routine uses XERROR
!     to print the following error message:  "SPL1D1 N LESS THAN 4.
!     RESULTS INCORRECT."
!
!     BOUNDARY CONDITIONS:
!
!        IOP(1) and IOP(2) are integers = 1,2,3,4 or 5.  By proper
!     combinations of these parameters, the user has a variety of
!     allowable end conditions.  The information associated with the end
!     conditions must be stored in W(1) and W(K) where K = (N-1)*IJ+1.
!
!        If IOP(1) = 1, then W(1) is the value of the second derivative
!     at X .  If IOP(1) = 2, then W(1) determines F''(X ) by the
!         1                                            1
!     relation F''(X ) = W(1) * F''(X ).  If IOP(1) = 3, then W(1) is
!                   1                2
!     the value of the first derivative at X .  IOP(1) = 4 allows
!                                           1
!     periodic boundary conditions, i.e., F'' = F'', F  = F .  If
!                                          1     N    1    N
!     IOP(1) = 5, the first derivative at X, is calculated by using a
!     differentiated four-point Lagrangian interpolation formula
!     evaluated at X .  For this case, W(1) need not be specified.
!                   1
!
!     The paragraph above describes the possibilities at the left
!     boundary.  The same holds at the right boundary, i.e., at X ,
!                                                                N
!     with the necessary boundary information stored at W(K). For
!     example, if IOP(2) = 2 then W(K) determines F''(X ) by the
!                                                      N
!     relation F''(X ) = W(K) * F''(X   ).
!                   N                N-1
!
!         Also, any combination of the above boundary conditions with
!     the exception of the periodic case may be used.  If the user
!     knows the second derivative at the left boundary and the first
!     derivative at the right boundary then IOP=(1) = 1, IOP(2) = 3,
!     W(1) = y'', W(K) = y'.  For LOGICAL results, the user must have
!             1           N
!     IOP(1) = 4 and IOP(2) = 4 for the periodic boundary conditions.
!
!     CAUTION:
!
!        In options 2 and 3, the input data in W(1) and W(K) is
!     destroyed by SPL1D1.  These boundary conditions must be reset if
!     SPL1D1 is to be used again.  For option 5, the user must have
!     N >= 4.
!
!
!     METHOD:
!
!        Given N knots a = X  < X  < X  ... X  = b on the interval
!                           1    2    3      N
!     [a,b] and prescribed function values at the knots, say
!     F(X ) = F , k = 1, 2, ..., N, cubics are joined between the points
!        k     k
!     (X ,F ) with the requirement that first and second derivatives
!       k  k
!     be continuous at these points.  This requirement yields N-2
!     equations in the N unknowns F'', F'', ..., F'' which when coupled
!                                  1    2         N
!     with two boundary conditions yields a system of equations solvable
!     for the quantities F''.  A complete discussion of the method may
!                         k
!     be found in the reference listed below.
!
!***REFERENCES  Walsh, Ahlberg, and Nilson, Journal of Mathematics and
!                 Mechanics, Vol. II, No. 2, 1962.
!***ROUTINES CALLED  XERROR
!***END PROLOGUE  SPL1D1
SAVE
INTEGER IOP(2), N, K, M, j1, j2, k2, k3, ij, i, k1, k4, j, l1, l2, l3, mk, ml
REAL(Kind=WP_Kind) con, don, e, c1, c2, a1, a2, a3, a4, a5, a6
REAL(Kind=WP_Kind) b1, b2, b3, b4, b5, b6, bob, bill, d1, d2
REAL(Kind=WP_Kind) X (N), F (N), W (N), A (N), B (N), C (N)
!***FIRST EXECUTABLE STATEMENT  SPL1D1
K = N - 1
A (2) = - (X (2) - X (1) ) / 6.d0
B (2) = (X (3) - X (1) ) / 3.d0
W (IJ + 1) = (F (2 * IJ + 1) - F (IJ + 1) ) / (X (3) - X (2) ) - (F (IJ + 1) - F (1) ) / (X (2) - X (1) )
IF(N==3)GOTO 4
   DO I = 3, K
      M = (I - 1) * IJ + 1
      J1 = M + IJ
      J2 = M - IJ
      CON = (X (I + 1) - X (I - 1) ) / 3.d0
      DON = (X (I) - X (I - 1) ) / 6.d0
      B (I) = CON - (DON**2) / B (I - 1)
      E = (F (J1) - F (M) ) / (X (I + 1) - X (I) ) - (F (M) - F (J2)) / (X (I) - X (I - 1) )
      W (M) = E- (DON * W (J2) ) / B (I - 1)
      A (I) = - (DON * A (I - 1) ) / B (I - 1)
   ENDDO
4  K1 = (N - 2) * IJ + 1
   C (N - 1) = - ( (X (N) - X (N - 1) ) / 6.d0) / B (N - 1)
   W (K1) = W (K1) / B (N - 1)
   A (N - 1) = A (N - 1) / B (N - 1)
   K2 = K - 1
IF(N==3)GOTO 8
   DO I = 2, K2
      J = N - I
      CON = (X (J + 1) - X (J) ) / 6.d0
      A (J) = (A (J) - CON * A (J + 1) ) / B (J)
      C (J) = - (CON * C (J + 1) ) / B (J)
      K3 = (J - 1) * IJ + 1
      M = K3 + IJ
      W (K3) = (W (K3) - CON * W (M) ) / B (J)
   ENDDO
8  K4 = (N - 1) * IJ + 1
IF(IOP (1)==5)GOTO 200
  201 C1 = W (1)
IF(IOP (2)==5)GOTO 202
  203 C2 = W (K4)
GOTO 205
  200 IF(N<4)GOTO 300
  302 A1 = X (1) - X (2)
A2 = X (1) - X (3)
A3 = X (1) - X (4)
A4 = X (2) - X (3)
A5 = X (2) - X (4)
A6 = X (3) - X (4)
W (1) = F (1) * (1.d0 / A1 + 1. / A2 + 1.d0 / A3) - A2 * A3 * F (IJ + &
 1) / (A1 * A4 * A5) + A1 * A3 * F (2 * IJ + 1) / (A2 * A4 * A6) &
 - A1 * A2 * F (3 * IJ + 1) / (A3 * A5 * A6)
GOTO 201
  202 IF(N<4)GOTO 300
  303 B1 = X (N) - X (N - 3)
B2 = X (N) - X (N - 2)
B3 = X (N) - X (N - 1)
B4 = X (N - 1) - X (N - 3)
B5 = X (N - 1) - X (N - 2)
B6 = X (N - 2) - X (N - 3)
L1 = K4 - IJ
L2 = L1 - IJ
L3 = L2 - IJ
W (K4) = - B2 * B3 * F (L3) / (B6 * B4 * B1) + B1 * B3 * F (L2) &
 / (B6 * B5 * B2) - B1 * B2 * F (L1) / (B4 * B5 * B3) + F (K4) &
 * (1.d0 / B1 + 1.d0 / B2 + 1.d0 / B3)
GOTO 203
  205 CONTINUE
I = 0
   50 CONTINUE
I = I + 1
IF(I>K) goto 100
!
M = (I - 1) * IJ + 1
GOTO 60
   70 IF(I==1)GOTO 50
   80 W (1) = W (1) - BOB * W (M)
W (K4) = W (K4) - BILL * W (M)
A (1) = A (1) - BOB * A (I)
A (N) = A (N) - BILL * A (I)
C (1) = C (1) - BOB * C (I)
C (N) = C (N) - BILL * C (I)
goto 50
!
   60 MK = IOP (1)
IF(MK==1)GOTO 62
IF(MK==2)GOTO 64
IF(MK==3)GOTO 66
IF(MK==4)GOTO 68
IF(MK==5)GOTO 66
   62 IF(I/=1)GOTO 71
   63 A (1) = - 1.d0
C (1) = 0.
GOTO 500
   71 BOB = 0.d0
GOTO 500
   64 IF(I/=1)GOTO 73
   76 A (1) = - 1.d0
C (1) = 0.
W (1) = 0.
GOTO 500
   73 IF(I>2)GOTO 82
   81 BOB = - C1
GOTO 500
   82 BOB = 0.d0
GOTO 500
   66 IF(I/=1)GOTO 83
   84 A (1) = - (X (2) - X (1) ) / 3.d0
C (1) = 0.d0
W (1) = - C1 + (F (IJ + 1) - F (1) ) / (X (2) - X (1) )
GOTO 500
   83 IF(I>2)GOTO 86
   85 BOB = (X (2) - X (1) ) / 6.d0
GOTO 500
   86 BOB = 0.d0
GOTO 500
   68 IF(I/=1)GOTO 87
   88 A (1) = - 1.d0
C (1) = 1.d0
W (1) = 0.d0
GOTO 500
   87 BOB = 0.d0
  500 ML = IOP (2)
IF(ML==1)GOTO 120
IF(ML==2)GOTO 130
IF(ML==3)GOTO 140
IF(ML==4)GOTO 150
IF(ML==5)GOTO 140
  120 IF(I/=1)GOTO 121
  122 A (N) = 0.d0
C (N) = - 1.d0
GOTO 70
  121 BILL = 0.d0
GOTO 70
  130 IF(I/=1)GOTO 131
  132 A (N) = 0.d0
C (N) = - 1.d0
W (K4) = 0.
GOTO 70
  131 IF(I/=K)GOTO 134
  133 BILL = - C2
GOTO 70
  134 BILL = 0.d0
GOTO 70
  140 IF(I/=1)GOTO 141
  142 A (N) = 0.d0
C (N) = (X (N - 1) - X (N) ) / 3.d0
W (K4) = C2 - (F (K4) - F (K1) ) / (X (N) - X (N - 1) )
GOTO 70
  141 IF(I/=K)GOTO 143
  144 BILL = (X (N) - X (N - 1) ) / 6.d0
GOTO 70
  143 BILL = 0.d0
GOTO 70
  150 IF(I/=1)GOTO 151
  152 A (N) = 0.
C (N) = (X (N - 1) + X (1) - X (N) - X (2) ) / 3.d0
W (K4) = (F (IJ + 1) - F (1) ) / (X (2) - X (1) ) - (F (K4) &
 - F (K1) ) / (X (N) - X (N - 1) )
GOTO 70
  151 IF(I/=2)GOTO 153
  154 BILL = (X (2) - X (1) ) / 6.d0
GOTO 70
  153 IF(I/=K)GOTO 155
  156 BILL = (X (N) - X (N - 1) ) / 6.d0
GOTO 70
  155 BILL = 0.
GOTO 70
  100 CON = A (1) * C (N) - C (1) * A (N)
D1 = - W (1)
D2 = - W (K4)
W (1) = (D1 * C (N) - C (1) * D2) / CON
W (K4) = (A (1) * D2 - D1 * A (N) ) / CON
DO I = 2, K
   M = (I - 1) * IJ + 1
   W (M) = W (M) + A (I) * W (1) + C (I) * W (K4)
ENDDO
GOTO 305
  300 CALL ErrWrt (' SPL1D1,  N<4 RESULTS INCORRECT')
  305 RETURN
ENDSUBROUTINE SPL1D1
