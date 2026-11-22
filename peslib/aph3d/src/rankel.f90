SUBROUTINE rankel (l, x, r, rp )
USE FileUnits_Module
!
IMPLICIT NONE
INTEGER l, i
REAL(Kind=WP_Kind) rx, f, rmax, fg, fgp, r, rp
REAL(Kind=WP_Kind) x, fgo, fgpo, fgn, fgpn, t
REAL(Kind=WP_Kind), PARAMETER::  zero=0.d0, one=1.d0 ,xmax=100.d0

rx = one / x
IF(x>xmax)THEN
f = exp ( -xmax )
ELSEIF(x<-xmax)THEN
f = exp ( xmax )
ELSE
f = exp ( -x )
ENDIF
!  IF(f==0.d0)THEN
!     WRITE(Out_Unit,*)'exponent too large: x=',x
!     STOP 'rankel'
!  ENDIF

fgo = f
fg = f*(one + rx)
fgpo = zero
fgp = -rx*rx*f

IF( l == 0 )THEN
    r = -fgo
    rp = -r
    RETURN
ELSEIF( l==1 )THEN
    r = -fg
    rp = fg-fgp
    RETURN
ELSE
    t = 3.d0
    i = 1
3    fgn = fg*t*rx + fgo
    fgpn = t*rx*(fgp - fg*rx) + fgpo
    i = i + 1
    t = t + 2.d0
    fgo = fg
    fgpo = fgp
    fg = fgn
    fgp = fgpn
    IF( i < l ) GOTO 3
    r = fg
    rp = fgp-fg
    RETURN
ENDIF

END
