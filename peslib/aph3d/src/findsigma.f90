SUBROUTINE findsigma(r, hstep, v0, v1, v2, r1, vl, r2, vr, vasy, jval, mu)
USE Numeric_Kinds_Module
USE Numbers_Module
USE FileUnits_OneDim_Module
IMPLICIT NONE
INTEGER jval
REAL(KIND=WP_Kind), PARAMETER:: tol=Ten*EPSILON(One)
REAL(KIND=WP_Kind) r, hstep, r1, vl, r2, vr, vasy, delta, v0, v1, v2, mu
LOGICAL, PARAMETER:: idbug=.false.

vl=vl-vasy
vr=vr-vasy
r=(r1+r2)/Two
delta=ABS(r1-r2)/r
CALL PotDiatomic(r, hstep, v0, v1, v2, jval, mu)
v0=v0-vasy
DO WHILE(delta.gt.tol)
  CALL Bisect(r, v0, r1, vl, r2, vr, delta)
  CALL PotDiatomic(r, hstep, v0, v1, v2, jval, mu)
  v0=v0-vasy
ENDDO
r=(r1+r2)/Two
CALL PotDiatomic(r, hstep, v0, v1, v2, jval, mu)
IF(idbug)THEN
   WRITE(Dbug_Unit,'(4(a,es23.15))')'sigma=', r,'  V(sigma)=', v0,'  V1(sigma)=', v1, '  V2(sigma)=', v2
ENDIF

RETURN
END SUBROUTINE FindSigma
