SUBROUTINE findextreme(i, nmax, rvals, hstep, vp, jval, mu, rext, v0, v1, v2, root)
USE Numeric_Kinds_Module
USE Numbers_Module
USE FileUnits_OneDim_Module
IMPLICIT NONE
LOGICAL root
INTEGER i, nmax, jval
REAL(KIND=WP_Kind) rext, hstep, mu, r1, r2, vl, vr, v0, v1, v2
REAL(KIND=WP_Kind) rvals(0:nmax), vp(0:nmax)

root=.false.
IF(vp(i-1)*vp(i).lt.Zero)THEN
   r1=rvals(i-1)
   vl=vp(i-1)
   r2=rvals(i)
   vr=vp(i)
   CALL findext(rext, hstep, v0, v1, v2, r1, vl, r2, vr, jval, mu)
   root=.true.
ENDIF

RETURN
END SUBROUTINE FindExtreme
