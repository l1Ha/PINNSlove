SUBROUTINE RhoCalc(Rhomin, DeltaRho1, DeltaRho2, IthRho, Rho, Rholast, Rhonext)
USE Numeric_Kinds_Module
USE FileUnits_Module
IMPLICIT NONE
LOGICAL, PARAMETER:: DEBUG=.FALSE.
INTEGER IthRho, N
REAL(Kind=WP_Kind) Rhomin, DeltaRho1, DeltaRho2, Rho, Rholast, Rhonext

N = IthRho - 1
Rho = (Rhomin + N*DeltaRho1)*(1.d0+DeltaRho2)**N

N = IthRho - 2
Rholast = (Rhomin + N*DeltaRho1)*(1.d0+DeltaRho2)**N

N = IthRho
Rhonext = (Rhomin + N*DeltaRho1)*(1.d0+DeltaRho2)**N

IF(DEBUG)THEN
   WRITE(Out_Unit,'("For IthRho = ",I5," Rholast = ",ES22.15, " Rho = ",ES22.15," Rhonext = ",ES22.15)')IthRho,Rholast,Rho,Rhonext
   WRITE(Out_Unit,'("Rhomin=",es23.15," DeltaRho1=",es23.15," DeltaRho2=",es23.15)')Rhomin,DeltaRho1,DeltaRho2
ENDIF

RETURN
ENDSUBROUTINE RhoCalc
