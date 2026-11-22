SUBROUTINE CenTerm(r, jval, mu, cent0, cent1, cent2, hstep)
USE Numeric_Kinds_Module
USE Numbers_Module
USE FileUnits_OneDim_Module
IMPLICIT NONE
INTEGER jval
REAL(KIND=WP_Kind) r, mu, cent0, cent1, cent2, hstep, rmincent
!REAL(KIND=WP_Kind), PARAMETER :: rmincent=1.d-10

rmincent=hstep
IF(jval.ne.0)THEN
   IF(r>zero)THEN
      cent0=jval*(jval+1)/(Two*mu*r**2)
      cent1=-Two*cent0/r
      cent2= Six*cent0/r**2
   ELSE
      cent0=jval*(jval+1)/(Two*mu*rmincent**2)
      cent1=-Two*cent0/rmincent
      cent2= Six*cent0/rmincent**2
   ENDIF
ELSE
   cent0=Zero
   cent1=Zero
   cent2=Zero
ENDIF

RETURN
END SUBROUTINE CenTerm
