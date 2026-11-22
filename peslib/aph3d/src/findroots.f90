SUBROUTINE findroots(i, nmax, rvals, hstep, v, vasy, jval, mu, sigma, v0, v1, v2, root)
USE Numeric_Kinds_Module
USE Numbers_Module
USE FileUnits_OneDim_Module
IMPLICIT NONE
LOGICAL root
INTEGER i, nmax, jval
REAL(KIND=WP_Kind) vasy, sigma, hstep, mu, r1, r2, vl, vr, v0, v1, v2
REAL(KIND=WP_Kind) rvals(0:nmax), v(0:nmax)

root=.false.
IF((v(i-1)-vasy)*(v(i)-vasy).lt.Zero)THEN
   r1=rvals(i-1)
   vl=v(i-1)
   r2=rvals(i)
   vr=v(i)
   CALL findsigma(sigma, hstep, v0, v1, v2, r1, vl, r2, vr, vasy, jval, mu)
   root=.true.
ENDIF

RETURN
END SUBROUTINE FindRoots
