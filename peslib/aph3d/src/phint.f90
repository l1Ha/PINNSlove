SUBROUTINE phint (phirms, ntheta, nchi, thetaval, chivals)
USE FileUnits_Module
USE Masses_Module
USE MassFactor2_Module
USE STST_Module
USE fbeta_Module
!
!     $RCSfile: phint.f,v $ $Revision: 1.14 $
!     $Date: 89/10/18 14:18:14 $
!     $State: Stable $
!

!     I N P U T    A R G U M E N T S
!
!     phirms
!     ntheta
!     nchi
!     thetaval
!     chivals

IMPLICIT NONE


!     I N T E G E R S
INTEGER i, ntheta, nchi, jchi, itheta

!     R E A L S
REAL(Kind=WP_Kind) thetaval, chivals, gtot, gtheta, gchi,  theta, chi, phirms, betasav


!     D I M E N S I O N S
DIMENSION phirms(ntheta,nchi), thetaval(ntheta), chivals(nchi)
DIMENSION betasav(narran)

!     E X T E R N A L S
EXTERNAL grms, angles


CALL angles (thetaval, chivals, ntheta, nchi)

DO 10 i=1,narran
   betasav(i)=betaw(i)
   betaw(i)=betafrst(i)
 10   CONTINUE
!-----------------------------------------------------------------------
DO 2 jchi = 1, nchi
   DO 1 itheta = 1, ntheta
      theta = thetaval(itheta)
      chi   = chivals(jchi)
      CALL grms (theta, chi, gtot, gtheta, gchi)
      phirms(itheta, jchi) = gtot
 1       CONTINUE
 2    CONTINUE

DO 20 i=1,narran
   betaw(i)=betasav(i)
 20   CONTINUE
RETURN
END
