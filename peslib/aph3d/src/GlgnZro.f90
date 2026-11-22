SUBROUTINE GlgnZro ( n, x )
USE Numeric_Kinds_Module
USE Numbers_Module
USE FileUnits_Module
!
!-------------------- GlgnZro -----------------------------------------
!
! GlgnZro -- zeroes of Guass Legendre polynomials
!
! === Latest Revision:    | 900524 at 15:38 | fixed up the convergence
!
!  input --
!
!   n -- order of the Guass Legendre polynomial whose zeroes are desired
!        (h0 has NONE, h1 has one zero, h25 has 25 zeroes)
!
!  output --
!
!     x -- array listing the non-negative zeroes of nth polynomial
!          x will be returned with 13 numbers in it IF n=25
!           note-- dimension statements limit n=200
!
!-----------------------------------------------------------------------
IMPLICIT NONE
SAVE
CHARACTER(LEN=21), PARAMETER:: ProcName='GlgnZro'
INTEGER, PARAMETER:: Maxlval=2001   ! PORT: enlarged from 201
LOGICAL   ODD, LAST
INTEGER i, ifpr, iodd, iroot, it, n, nnr, np, npr
REAL(Kind=WP_Kind) xg(Maxlval), hp(Maxlval), x(*)
REAL(Kind=WP_Kind) deriv, cor, theta, xz, xzn
!
!------------------- execution begins here ----------------------------
!
IF( n < 1 ) RETURN
!
x(1) = 0d0
IF( n == 1 ) RETURN
!
IF(n+1>Maxlval)THEN
   WRITE(Out_Unit,*)'Error in GlgnZro :',n+1, Maxlval
   WRITE(Msg_Unit,*)'Error in GlgnZro :',n+1, Maxlval
   STOP 'GlgnZro'
ENDIF

npr = n/2
np = n + 1
iodd = mod(n,2)
ODD = iodd == 1
nnr = iodd + npr
ifpr = iodd + 1
IF( ODD ) x(1) = 0d0
!
DO iroot = ifpr, nnr
   theta = pi*(4.d0*(npr+iroot+1)-5.d0)/(4.0d0*n+2.0d0)
   theta = theta + 1.d0/(8.d0*n**2*tan(theta))
   xg(iroot) = -cos(theta)
ENDDO
!
DO i = ifpr, nnr
   xz = xg(i)
   CALL Pn ( xz, n, hp)
   
   it = 0
   LAST = .false.
!
 11      CONTINUE
      it = it + 1
      deriv = (n*(-xz*hp(np)+hp(n)))/(1.d0-xz*xz)
      cor = hp(np)/deriv
      xzn = xz - cor
!
      IF( LAST ) GOTO 12
      !LAST = 1d0 + .5d0*ABS(cor/xzn) == 1d0
      LAST=ABS(0.5d0*abs(cor/xzn))<1.d-12
      IF( it > 30 ) GOTO 13
!
      CALL Pn ( xzn, n, hp)
      xz = xzn
   GOTO 11
!
 13      CONTINUE
   WRITE(Out_Unit, 1000 ) it, iroot, n, xzn, xz
   WRITE(Msg_Unit, 1000 ) it, iroot, n, xzn, xz
   STOP 'GlgnZro'
!
 12      CONTINUE
   x(i) = xzn
ENDDO
!
 1000 FORMAT( ' poor convergence after ', i2, ' iterations *** GlgnZro -- root(order) = '  &
             , i3,"(",i3,")",/,' best zero = ', e20.10, ', previous guess = ', e20.10 )
!
RETURN
ENDSUBROUTINE GlgnZro
