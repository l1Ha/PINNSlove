SUBROUTINE devcsf (na, hamil, nb, eig, vec, nc)
USE Numeric_Kinds_Module
USE FileUnits_Module
IMPLICIT NONE
INTEGER :: na, nb, nc, info, lwork
REAL(Kind=WP_Kind) :: hamil (na, na), vec (na, na), eig (na)
REAL(Kind=WP_Kind) :: work (200, 200)

lwork = 200 * 200
CALL dsyev ('V', 'U', na, hamil, na, eig, work, lwork, info)
vec = hamil

RETURN
ENDSUBROUTINE devcsf
