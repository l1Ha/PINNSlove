FUNCTION Indxr_New(lam, lammin, nrho, ntheta, nchi, irho, itheta, ichi)
IMPLICIT NONE

INTEGER :: Indxr_New, nrho, ntheta, nchi, irho, itheta, ichi, lam, lammin
!c      Indxr_New=ichi+(itheta-1)*nchi+(irho-1)*nchi*ntheta
!c     >     +(lam-lammin)*nchi*ntheta*nrho

Indxr_New = irho+(itheta-1)*nrho+(ichi-1)*nrho*ntheta+(lam-lammin)*nchi*ntheta*nrho
RETURN
END FUNCTION Indxr_New
