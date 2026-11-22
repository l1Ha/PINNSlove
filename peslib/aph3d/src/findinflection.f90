SUBROUTINE findinflection(i, nmax, rvals, hstep, vpp, jval, mu, rext, v0, v1, v2, root)
USE Numeric_Kinds_Module
USE Numbers_Module
USE FileUnits_OneDim_Module
IMPLICIT NONE
LOGICAL root
INTEGER i, nmax, jval
REAL(KIND=WP_Kind) rext, hstep, mu, r1, r2, vl, vr, v0, v1, v2
REAL(KIND=WP_Kind) rvals(0:nmax), vpp(0:nmax)

root=.false.
IF(vpp(i-1)*vpp(i).lt.Zero)THEN
   r1=rvals(i-1)
   vl=vpp(i-1)
   r2=rvals(i)
   vr=vpp(i)
   CALL findinf(rext, hstep, v0, v1, v2, r1, vl, r2, vr, jval, mu)
   root=.true.
ENDIF

RETURN
END SUBROUTINE FindInflection
