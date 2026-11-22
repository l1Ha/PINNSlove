SUBROUTINE legndr(jmax, nbtheta, pj, xpt)
USE Numeric_Kinds_Module
USE Numbers_Module

! evaluates associated legendre polynomials at the quadrature
! points.

IMPLICIT NONE
CHARACTER(LEN=21), PARAMETER:: ProcName='legndr  '
INTEGER nbtheta, ibtheta, jmax, j, lambda
REAL(Kind=WP_Kind) xpt(nbtheta), pj(nbtheta,0:jmax,0:jmax), plm_fun
REAL(Kind=WP_Kind) sqrt2


DO ibtheta=1,nbtheta
   DO j=0,jmax
      DO lambda=0,j
         pj(ibtheta,j,lambda)=plm_fun(j, lambda, xpt(ibtheta))
      ENDDO
   ENDDO
ENDDO

! On return array pj contains associate legendre polynomial
! values at the quadrature points.

RETURN
END
