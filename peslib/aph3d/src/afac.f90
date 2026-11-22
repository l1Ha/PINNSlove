FUNCTION afac (jtot, lam, lamp, parity)
USE Numeric_Kinds_Module
USE FileUnits_Module
IMPLICIT NONE
INTEGER :: jtot, lam, lamp, parity

REAL(Kind=WP_Kind) :: afac, lambda_plus, lambda_minus, delta
!---------------------------------------------------------------------
! Factors for Asymmetric Top coupling term. See APH theory paper
! Eq. (176).
!---------------------------------------------------------------------

afac = 0.5d0 / sqrt ( (1.d0 + delta (lam, 0) ) * (1.d0 + delta ( &
 lamp, 0) ) ) * (lambda_plus (jtot, lam) * lambda_plus (jtot, lam + &
 1) * delta (lamp, lam + 2) + lambda_minus (jtot, lam) * &
 lambda_minus (jtot, lam - 1) * delta (lamp, lam - 2) + ( - 1) ** ( &
 jtot + lam + parity) * lambda_minus (jtot, lam) * lambda_minus ( &
 jtot, lam - 1) * delta (lamp, 2 - lam) )
RETURN
ENDFUNCTION afac
