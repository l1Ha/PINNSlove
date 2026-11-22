SUBROUTINE lgnd11 (na, np, w, x)
USE Numeric_Kinds_Module
USE Numbers_Module
USE FileUnits_Module
INTEGER na, np, i
REAL(Kind=WP_Kind) w(na), x(na)
!  calculates points and weights needed in Gauss_Legendre quadrature on
!  na is number of angles used; np is number of points in quadrature
!  for symmetric systems one can use np=2*na; otherwise np=na.
!  requires SUBROUTINE lgndrx
!  DO 10 i=1, na
!  CALL lgndrx (np, i, w(i), x(i))
!  convert from (0,1) points and weights of lgndrx to those for (-1,1).
!  x(i)=2.*x(i)-1.
!  w(i)=2.*w(i)
!  10 CONTINUE
CALL glegen(np, x, w, -One, One)
WRITE(Out_Unit,20) (x(i), i=1, na)
WRITE(Out_Unit,30) (w(i), i=1, na)
RETURN
!
 20   FORMAT(/, 'x(i)=cos=', 11f10.7)
 30   FORMAT(/, 'weight(i)=', 11f10.7)
ENDSUBROUTINE Lgnd11
