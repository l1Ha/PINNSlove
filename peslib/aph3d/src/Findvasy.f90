SUBROUTINE Findvasy(vasy, ri, hstep)
USE Numeric_Kinds_Module, ONLY: WP_Kind, IW_Kind
USE Numbers_Module,       ONLY: WP_Epsilon, Two
USE OneDim_Parms_Module,  ONLY: IARRAN
USE FileUnits_Common_Module
USE FileUnits_OneDim_Module
IMPLICIT NONE
INTEGER(KIND=IW_Kind) i
REAL(KIND=WP_Kind)    vasy, ri, hstep, v0, v1, v2, mu, r1, r2, vl, vr
REAL(KIND=WP_Kind), PARAMETER:: tol=1.d-7
REAL(KIND=WP_Kind), PARAMETER :: Power10=1.E10_WP_Kind

CALL TimeInfoSub('Findvasy','Begin Findvasy')
CALL PotDiatomic(Power10, hstep, vasy, v1, v2, 0, mu)
WRITE(Out_Unit,*)"vasy=",vasy
DO i=0,10
   ri=Two**i
   CALL PotDiatomic(ri, hstep, v0, v1, v2, 0, mu)
   WRITE(Out_Unit,'(A,es23.15,A,es23.15)')"ri= ", ri, "  v0-vasy=", v0-vasy
   WRITE(0,'(A,es23.15,A,es23.15)')"ri= ", ri, "  v0-vasy=", v0-vasy
   IF(ABS(v0-vasy).le.tol)THEN
      r1=Two**(i-1)
      CALL PotDiatomic(r1, hstep, vl, v1, v2, 0, mu)
      r2=Two**i
      CALL PotDiatomic(r2, hstep, vr, v1, v2, 0, mu)
      CALL findri(ri, hstep, v0, v1, v2, r1, vl, r2, vr, vasy, 0, mu, tol)
      EXIT
   ENDIF
ENDDO
CALL TimeInfoSub('Findvasy','  End Findvasy')

RETURN
END SUBROUTINE Findvasy
