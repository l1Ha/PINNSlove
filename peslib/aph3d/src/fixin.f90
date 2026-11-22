SUBROUTINE fixin (idt, ntheta, nchi, thetaval, chivals)
USE Numeric_Kinds_Module
USE FileUnits_Module
USE FEM_Module
!
!     $RCSfile: fixin.f,v $ $Revision: 1.14 $
!     $Date: 89/10/18 14:18:03 $
!     $State: Stable $
!

!     I N P U T    A R G U M E N T S
!
!     idt      Identifies each node as free, fixed, constrained, unique, nonexistant
!     ntheta
!     nchi
!     thetaval
!     chivals


IMPLICIT NONE
INTEGER idt, ntheta, nchi, jchi, itheta
REAL(Kind=WP_Kind) thetaval, chivals, tolgtot, gtot, gtheta, gchi, theta, chi
DIMENSION idt(ntheta,nchi), thetaval(ntheta), chivals(nchi)
EXTERNAL grms, angles
DATA tolgtot/1.d-100/

CALL angles (thetaval, chivals, ntheta, nchi)

DO jchi = 1, nchi
   DO itheta = 1, ntheta
      IF(idt(itheta,jchi)==free)THEN
         theta = thetaval(itheta)
         chi   = chivals(jchi)
         CALL grms(theta, chi, gtot, gtheta, gchi)
         IF(gtot<tolgtot)THEN
            idt(itheta,jchi) = fixed
         ENDIF
      ENDIF
   ENDDO
ENDDO

RETURN
END
