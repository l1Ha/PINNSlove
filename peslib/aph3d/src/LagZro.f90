SUBROUTINE LagZro ( n, x )
USE Numeric_Kinds_Module
USE Numbers_Module
!
!-------------------- LagZro -----------------------------------------
!
! LagZro -- zeroes of Laguere polynomials
!
! === Latest Revision:    | 900524 at 15:38 | fixed up the convergence
!
!  input --
!
!     n -- order of the Laguere polynomial whose zeroes are desired
!          (h0 has NONE, h1 has one zero, h25 has 25 zeroes)
!
!  output --
!
!     x -- array listing the non-negative zeroes of nth polynomial
!          x will be returned with 13 numbers in it IF n=25
!           note-- dimension statements limit n=200
!
!---------------------------------------------------------------------
!
USE FileUnits_Module

IMPLICIT NONE
!
!
LOGICAL   LAST
INTEGER   n, np, it, i, ibisect
!
REAL(Kind=WP_Kind) xg(300), hp(300), x(n), lower(300), upper(300)
REAL(Kind=WP_Kind) beta, beta8, beta8sq, bj, max, delta, xl, xu, xz, cor, abs, xzn, hl, hu
!
!------------------ execution begins here ----------------------------
!
IF( n < 1 ) RETURN
!
x(1) = 1.0d0
IF( n == 1 ) RETURN
!
np = n + 1
!
DO 1 i = 1, n
   beta = (i-0.25d0)*pi
   beta8 = 8.d0*beta
   beta8sq = beta8**2
   bj = beta + (1.d0 + (-124.d0/3.d0 + (120928.d0/15.d0 - 401743168.d0/105.d0/beta8sq)/beta8sq)/beta8sq)/ beta8
   xg(i) = bj**2*(1.d0+(-2.d0+bj**2)/(48.d0*(n+0.5d0)**2))/ (4.d0*(n+0.5d0))
   lower(i) = bj**2/(4.d0*(n+0.5d0))
   upper(i) = (i+0.5d0)/(n+0.5d0)*(2.d0*(i+0.5d0)+sqrt(4.d0*(i+0.5d0)**2+0.25d0))
 1    CONTINUE
!
DO 10 i = 1, n
!
! Make lower bound greater than previous root.
!
   IF(i>2)THEN
      IF(lower(i)<x(i-1)+0.5d0*(x(i-1)-x(i-2))/2)THEN
         lower(i) = x(i-1)+0.5d0*(x(i-1)-x(i-2))/2
         IF(xg(i)<lower(i))xg(i)=(lower(i)+upper(i))/2
      ENDIF
   ENDIF
   xz = xg(i)
   CALL Laguer_fun ( xz, n, hp, np )
   it = 0
   LAST = .false.
!
 11      CONTINUE
      it = it + 1
      cor = xz/n*hp(np)/(hp(np)-hp(n))

      IF(it==1.or.ABS(cor)>0.5d0*(upper(i)-lower(i)).or.xz-cor<lower(i).or.xz-cor>upper(i))THEN
         delta = (upper(i)-lower(i))/10
         DO 30 ibisect = 1, 10
            xl = lower(i) + delta*(ibisect-1)
            CALL Laguer_fun (xl, n, hp, np)
            hl = hp(np)
            xu = lower(i) + delta*ibisect
            CALL Laguer_fun (xu, n, hp, np)
            hu = hp(np)
            IF(hl*hu<0.d0)THEN
               xzn = (xl*hu-xu*hl)/(hu-hl)
               cor = xz - xzn
               lower(i) = xl
               upper(i) = xu
               GOTO 31
            ENDIF
   30          CONTINUE
   31          CONTINUE
      ENDIF

      xzn = xz - cor

!
      IF( LAST ) GOTO 12
      LAST = 1d0 + .01d0*ABS(cor/max(1.d0,xzn)) == 1d0
      IF( it > 40 ) GOTO 13
!
      CALL Laguer_fun ( xzn, n, hp, np )
      xz = xzn
   GOTO 11
!
 13      CONTINUE
   WRITE(Msg_Unit, 1000 ) it, i, n, xzn, xz
   STOP 'lagzro'
!
 12      CONTINUE
   x(i) = xzn
 10   CONTINUE
!
 1000 FORMAT( ' poor convergence after ', i2, ' iterations ***' ' LagZro -- root(order) = ', i3,"(",i3,")",/,   &
               ' best zero = ', e20.10, ', previous guess = ', e20.10 )
!
RETURN
END
