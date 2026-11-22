SUBROUTINE phot_read (n, omega)
!---------------------------------------------------------------------
! reads file fph_unit - photodissociation inhomogeneity omega
!---------------------------------------------------------------------
USE FileUnits_MODULE
IMPLICIT NONE
!#include <Inc.unitno>
INTEGER n
REAL(Kind=WP_Kind) omega(n)

READ(fph_unit) omega

RETURN
END
