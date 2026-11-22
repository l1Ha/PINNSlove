SUBROUTINE findri(ri, hstep, v0, v1, v2, r1, vl, r2, vr, vasy, jval, mu, tol)
USE Numeric_Kinds_Module
USE Numbers_Module
USE FileUnits_Common_Module
USE FileUnits_OneDim_Module
IMPLICIT NONE
INTEGER jval
REAL(KIND=WP_Kind) :: tol
REAL(KIND=WP_Kind) ri, hstep, r1, vl, r2, vr, vasy, delta, v0, v1, v2, mu

ri=(r1+r2)/Two
delta=ABS(r1-r2)/ri
CALL PotDiatomic(ri, hstep, v0, v1, v2, jval, mu)
DO WHILE(delta.gt.tol)
   IF(ABS(v0-vasy).gt.tol)THEN
      r1=ri
      vl=v0
   ELSE
      r2=ri
      vr=v0
   ENDIF
   ri=(r1+r2)/Two
   delta=ABS(r1-r2)/ri
   CALL PotDiatomic(ri, hstep, v0, v1, v2, jval, mu)
   v0=v0
ENDDO
ri=(r1+r2)/Two
CALL PotDiatomic(ri, hstep, v0, v1, v2, jval, mu)

WRITE(Out_Unit,*)
WRITE(Out_Unit,'(A,es23.15)')"Found Asymptotic distance where V(r)<", tol
WRITE(Out_Unit,'(4(a,es23.15))')'ri=', ri,'  V(ri)=', v0,'  V1(ri)=', v1, '  V2(ri)=', v2

RETURN
END SUBROUTINE FindRi
