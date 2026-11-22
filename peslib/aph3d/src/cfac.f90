FUNCTION cfac (jtot, lam, lamp, parity)
USE Numeric_Kinds_Module
USE FileUnits_Module
IMPLICIT NONE
INTEGER :: jtot, lam, lamp, parity

REAL(Kind=WP_Kind) :: cfac, lambda_plus, lambda_minus, delta
!------------------------------------------------------------------
! Factors for the Coriolis coupling term. See Eq. (174) in the
! APH theory paper.
!------------------------------------------------------------------

cfac = 1.d0 / sqrt ( (1.d0 + delta (lam, 0) ) * (1.d0 + delta ( &
 lamp, 0) ) ) * (lambda_plus (jtot, lam) * delta (lamp, lam + 1) &
 - lambda_minus (jtot, lam) * delta (lamp, lam - 1) + lambda_minus &
 (jtot, lam) * ( - 1) ** (jtot + lam + parity) * delta (lamp, 1 - &
 lam) )
RETURN
END FUNCTION cfac
