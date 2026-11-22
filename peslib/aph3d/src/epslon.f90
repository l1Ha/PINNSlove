FUNCTION epslon (x)
USE Numeric_Kinds_Module
USE FileUnits_Module
IMPLICIT NONE
!
!  estimate unit roundoff in quantities of size x.
!
REAL(Kind=WP_Kind) x, a, b, c, eps, abs, epslon
INTRINSIC abs
!
!  this program should function properly on all systems
!  satisfying the following two assumptions,
!  1.  the base used in representing floating point
!  numbers is not a power of three.
!  2.  the quantity  a  in statement 10 is represented to
!  the accuracy used in floating point variables
!  that are stored in memory.
!  the statement number 10 and the GOTO 10 are intended to
!  force optimizing compilers to generate code satisfying
!  assumption 2.
!  under these assumptions, it should be true that,
!  a  is not exactly equal to four-thirds,
!  b  has a zero for its last bit or digit,
!  c  is not exactly equal to one,
!  eps  measures the separation of 1.0 from
!  the next larger floating point number.
!  the developers of eispack would appreciate being informed
!  about any systems where these assumptions DO not hold.
!
!  this version dated 4/6/83.
!
a = 4.0d0/3.0d0
 10   b = a - 1.0d0
c = b + b + b
eps = ABS(c-1.0d0)
IF(eps == 0.0d0) GOTO 10
epslon = eps*ABS(x)

RETURN
END FUNCTION epslon
