SUBROUTINE HrmZro ( n, x )
!
!---------------------- HrmZro -----------------------------------------
!
! HrmZro -- zeroes of Hermite polynomials
!
! === Latest Revision:    | 900524 at 15:38 | fixed up the convergence
!
!  input --
!
!     n -- order of the Hermite polynomial whose zeroes are desired
!          (h0 has NONE, h1 has one zero, h25 has 25 zeroes)
!
!  output --
!
!     x -- array listing the non-negative zeroes of nth polynomial
!          x will be returned with 13 numbers in it IF n=25
!           note-- dimension statements limit n=200
!
!-----------------------------------------------------------------------
!
USE Numeric_Kinds_Module
USE FileUnits_Module
USE Numbers_Module
INTEGER i, ic, ifpr, iodd, iroot, it, n, nnr, np, npr
REAL(Kind=WP_Kind) cor, xav, xfp, xs, xz, xzn, tn, tnr
REAL(Kind=WP_Kind) xg(200), hp(200), x(n)

LOGICAL   ODD, LAST
!
CHARACTER(LEN=2)  IGlobal
!
!---------------------xecution begins here ----------------------------
!

IF( n < 1 ) RETURN
!
x(1) = 0d0
IF( n == 1 ) RETURN
!
npr = n/2
np = n + 1
tn = 2d0*n + 1d0
tnr = .5d0/n
iodd = mod(n,2)
ODD = iodd == 1
nnr = iodd + npr
ifpr = iodd + 1
IF( ODD ) x(1) = 0d0
!
xfp = pi/sqrt(tn)
IF( .NOT.ODD ) xfp = 0.5d0*xfp
IF(      ODD ) xfp = pi/sqrt(tn-.25d0*xfp**2)
xg(1) = xfp
!
IF( npr < 2 ) GOTO 2
DO i = 2, npr
   xs    = pi/sqrt(tn - xg(i-1)**2)
   xav   = 0.5d0*xs + xg(i-1)
   xg(i) = xg(i-1) + pi/sqrt(tn - xav**2)
ENDDO

2 CONTINUE
!
ic = 0
iroot = npr + iodd
!
DO i = ifpr, nnr
   iroot = iroot + 1
   ic = ic + 1
   xz = xg(ic)
   CALL Hermit_fun ( xz, n, hp, np )
   it = 0
   LAST = .false.
!
11 CONTINUE
      it = it + 1
      cor = hp(np)*tnr/hp(n)
      xzn = xz - cor
!
      IF( LAST ) GOTO 12
      LAST = 1d0 + .5d0*ABS(cor/xzn) == 1d0
      IF( it > 30 ) GOTO 13
!
      CALL Hermit_fun ( xzn, n, hp, np )
      xz = xzn
   GOTO 11
!
13 CONTINUE
   WRITE(Msg_Unit, 1000 ) it, iroot, n, xzn, xz
!
12 CONTINUE
   x(i) = xzn
ENDDO
!
 1000 FORMAT( ' poor convergence after ', i2, ' iterations *** HrmZro -- root(order) = '  &
                  , i3,"(",i3,")",/, ' best zero = ', e20.10, ', previous guess = ', e20.10 )
!
RETURN
ENDSUBROUTINE HrmZro
