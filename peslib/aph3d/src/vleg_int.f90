REAL(Kind=WP_Kind) FUNCTION Vleg_Int(temp, vprime)
USE Numeric_Kinds_Module
USE numbers_Module
USE int1d_Module
IMPLICIT NONE
CHARACTER(LEN=21), PARAMETER:: ProcName='vleg_int'

! integration or quadrature routine. Integrates the PES over
! big-theta.

INTEGER k
REAL(Kind=WP_Kind) temp(nbtheta), vprime(nbtheta)

! perform Gauss_Legendre quadrature.

vleg_int=zero
DO k=1,nbtheta
   vleg_int=vleg_int+temp(k)*vprime(k)
ENDDO

RETURN
END FUNCTION Vleg_Int
