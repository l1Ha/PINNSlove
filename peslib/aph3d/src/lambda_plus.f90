FUNCTION lambda_plus (jtot, lambda)
USE Numeric_Kinds_Module
IMPLICIT NONE
INTEGER :: jtot, lambda

REAL(Kind=WP_Kind) :: lambda_plus
!-----------------------------------------------------------------------
! See APH theory paper Eq. (158)
!-----------------------------------------------------------------------
lambda_plus = (jtot + lambda + 1) * (jtot - lambda)

lambda_plus = sqrt (lambda_plus)
RETURN
END FUNCTION lambda_plus
