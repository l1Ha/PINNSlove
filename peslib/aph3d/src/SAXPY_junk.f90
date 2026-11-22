SUBROUTINE saxpy_junk(N,SA,SX,INCX,SY,INCY)
USE FileUnits_Module
!   BEGIN PROLOGUE  saxpy_junk
!   DATE WRITTEN   791001   (YYMMDD)
!   REVISION DATE  861211   (YYMMDD)
!   CATEGORY NO.  D1A7
!   KEYWORDS  LIBRARY=SLATEC(BLAS),
!             TYPE=SINGLE PRECISION(saxpy_junk-S DAXPY-D CAXPY-C),
!             LINEAR ALGEBRA,TRIAD,VECTOR
!   AUTHOR  LAWSON, C. L., (JPL)
!           HANSON, R. J., (SNLA)
!           KINCAID, D. R., (U. OF TEXAS)
!           KROGH, F. T., (JPL)
!   PURPOSE  S.P. computation y = a*x + y
!   DESCRIPTION
!
!                B L A S  Subprogram
!    Description of Parameters
!
!     --Input--
!        N  number of elements in input vector(s)
!       SA  single precision scalar multiplier
!       SX  single precision vector with N elements
!     INCX  storage spacing between elements of SX
!       SY  single precision vector with N elements
!     INCY  storage spacing between elements of SY
!
!     --Output--
!       SY  single precision result (unchanged IF N <= 0)
!
!     OverWRITE single precision SY with single precision SA*SX +SY.
!     For I = 0 to N-1, replace  SY(LY+I*INCY) with SA*SX(LX+I*INCX) +
!       SY(LY+I*INCY), where LX = 1 IF INCX >= 0, ELSE LX = (-INCX)*N
!       and LY is defined in a similar way using INCY.
!   REFERENCES  LAWSON C.L., HANSON R.J., KINCAID D.R., KROGH F.T.,
!                 *BASIC LINEAR ALGEBRA SUBPROGRAMS FOR FORTRAN USAGE*,
!                 ALGORITHM NO. 539, TRANSACTIONS ON MATHEMATICAL
!                 SOFTWARE, VOLUME 5, NUMBER 3, SEPTEMBER 1979, 308-323
!   ROUTINES CALLED  (NONE)
!   END PROLOGUE  saxpy_junk
!
INTEGER N, INCX, INCY, IX, IY, I, M, MP1, NS
REAL(Kind=WP_Kind) SX(1),SY(1),SA
!   FIRST EXECUTABLE STATEMENT  saxpy_junk
IF(N<=0.OR.SA==0.D0) RETURN
IF(INCX==INCY)THEN
 !IF(INCX-1) 5,20,60
 IF(INCX<1)GOTO 5
 IF(INCX==1)GOTO 20
 IF(INCX>1)GOTO 60
ENDIF
    5 CONTINUE
!
!        CODE FOR NONEQUAL OR NONPOSITIVE INCREMENTS.
!
IX = 1
IY = 1
IF(INCX<0)IX = (-N+1)*INCX + 1
IF(INCY<0)IY = (-N+1)*INCY + 1
DO 10 I = 1,N
  SY(IY) = SY(IY) + SA*SX(IX)
  IX = IX + INCX
  IY = IY + INCY
   10 CONTINUE
RETURN
!
!        CODE FOR BOTH INCREMENTS EQUAL TO 1
!
!
!        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 4.
!
   20 M = MOD(N,4)
IF( M == 0 ) GO TO 40
DO 30 I = 1,M
  SY(I) = SY(I) + SA*SX(I)
   30 CONTINUE
IF( N < 4 ) RETURN
   40 MP1 = M + 1
DO 50 I = MP1,N,4
  SY(I) = SY(I) + SA*SX(I)
  SY(I + 1) = SY(I + 1) + SA*SX(I + 1)
  SY(I + 2) = SY(I + 2) + SA*SX(I + 2)
  SY(I + 3) = SY(I + 3) + SA*SX(I + 3)
   50 CONTINUE
RETURN
!
!        CODE FOR EQUAL, POSITIVE, NONUNIT INCREMENTS.
!
   60 CONTINUE
NS = N*INCX
    DO 70 I=1,NS,INCX
    SY(I) = SA*SX(I) + SY(I)
   70     CONTINUE
RETURN
END
