Subroutine Basis_Jacobi(z, az, bz, cz)
USE fileunits_Module
USE numeric_kinds_Module
USE gaussb_Module
USE Masses_Module
USE Arrch_Module
USE DAS_Module
USE Numbers_Module
IMPLICIT NONE
REAL(Kind=WP_Kind)  z, az, bz, cz, g1, g2, g3, sm, kf, beta

kf   = weau(karran)/wexeau(karran)
beta = sqrt(two*usys*wexeau(karran))
sm   = re(karran)+one/beta*LOG(kf/(kf-one))
g2   = (beta)**2*kf/two*exp(-beta*(sm-re(karran)))
g3   = -(beta)**3*kf/two*exp(-beta*(sm-re(karran)))
az   = sqrt(g2)+g3*sm/six/sqrt(g2)
bz   = -g3*sm**3/six/sqrt(g2)
cz   = bz/sm-az*sm

RETURN
ENDSUBROUTINE Basis_Jacobi

