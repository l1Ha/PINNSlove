FUNCTION lambda_minus (Jtot, Lambda)
!-----------------------------------------------------------------
!   This Function Routine returns the value of Lambda_minus as
!   given in Equation (158) of the APH Theory paper. 
!
!   References
!   R. T Pack and G. A. Parker, J. Chem. Phys. 87 3988-3921 (1987).
! 
!   Required Input <=====
!   Jtot     Total angular momentum of the system
!   Lambda   Projection of the total angular momentum on the 
!            body-fixed z-axis.
!
!   On Return =====>
!   lambda_minus  The value of lambda_minus
!-----------------------------------------------------------------
USE Numeric_Kinds_Module
USE FileUnits_Module
IMPLICIT NONE
INTEGER, INTENT(IN) :: Jtot, Lambda

REAL(Kind=WP_Kind) :: lambda_minus
!-----------------------------------------------------------------
! See APH theory paper Eq. (158).
!-----------------------------------------------------------------
lambda_minus = (Jtot - Lambda + 1) * (Jtot + Lambda)

lambda_minus = SQRT(Lambda_minus)
RETURN
ENDFUNCTION lambda_minus
