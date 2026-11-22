SUBROUTINE Vleg_Mat(pj,vprime,weight,vleg,jmax,jtot,j_min,symmetry)
USE Numeric_Kinds_Module
USE int1d_Module
IMPLICIT NONE
LOGICAL symmetry
CHARACTER(LEN=21), PARAMETER:: ProcName='vleg_mat'

! integrate the PES(times associated legendre polynomials)
! over big-theta. 

INTEGER jf, jn, jmax, lambda, iltheta, jtot, k, j_min, deltaj
REAL(Kind=WP_Kind) pj(nbtheta,0:jmax,0:jmax)
REAL(Kind=WP_Kind) weight(nbtheta), temp(nbtheta) 
REAL(Kind=WP_Kind) vprime(nbtheta,nltheta+1), vleg_int
REAL(Kind=WP_Kind) vleg(nltheta+1,0:jmax,0:jmax,0:jmax)

! This corresponds to rtp98/4/30-5 Eq. 26. Big-theta integration
! only and hence the elements vleg are functions of little-theta.
deltaj=1
IF(symmetry)deltaj=2
DO jf=j_min,jmax,deltaj
   DO jn=j_min,jf,deltaj
      DO lambda=0,MIN(jf,jn,jtot)
         DO k=1,nbtheta
            temp(k) = pj(k,jf,lambda)*pj(k,jn,lambda)*weight(k)
         ENDDO
         DO iltheta=1,nltheta+1
               vleg(iltheta,jf,jn,lambda) = vleg_int(temp,vprime(1,iltheta))
               vleg(iltheta,jn,jf,lambda) = vleg(iltheta,jf,jn,lambda)
         ENDDO
      ENDDO
   ENDDO
ENDDO

! results are stored in vleg.

RETURN
ENDSUBROUTINE Vleg_Mat
