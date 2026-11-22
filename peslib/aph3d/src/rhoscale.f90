SUBROUTINE rhoscale(rho, wsub, wtemp, n1, n2, rhoval, asym)
USE Numeric_Kinds_Module
!
!-----------------------------------------------------------------------
! asym = LOGICAL variable .true. if wsub = asymmetric matrix
!-----------------------------------------------------------------------
IMPLICIT NONE
LOGICAL asym
INTEGER j, i, n1, n2, ij
REAL(Kind=WP_Kind) wsub(n1,n2), wtemp(n1*n2), rhoval, rho, rscale

ij=0
!-----------------------------------------------------------------------
! Symmetric Matrix
!-----------------------------------------------------------------------
IF(.NOT.asym)THEN
   rscale=(rhoval/rho)**2
   DO j=1,n1
      DO i = 1,j
         ij=ij+1
         wtemp(ij) = wtemp(ij)+rscale*wsub(i,j)
      ENDDO
   ENDDO
ELSE
!-----------------------------------------------------------------------
! Asymmetric Matrix
!-----------------------------------------------------------------------
   DO j=1,n2
      DO i = 1,n1
         ij=ij+1
         wtemp(ij) = wsub(i,j)/(rho*rho)
      ENDDO
   ENDDO
ENDIF

RETURN
END
