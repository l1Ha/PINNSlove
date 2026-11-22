SUBROUTINE angles (thetaval, chivals, ntheta, nchi)
USE FileUnits_Module
USE main_Module
USE delang_Module
USE mdata_Module
USE Numbers_Module
!


!  P U R P O S E   O F    S U B R O U T I N E
!-----------------------------------------------------------------------
!  This routine generates a list of theta and chi values that are used
!  in the finite elements.
!-----------------------------------------------------------------------

!  I N P U T    A R G U M E N T S

!  thetaval
!  chivals
!  ntheta
!  nchi


IMPLICIT NONE
INTEGER ntheta, nchi, k, nchidiv, jchimin, jchimax, kk, itheta, jchi, kthdiv, kkm1
REAL(Kind=WP_Kind)  thetaval, chivals, chistart, theta, chi
DIMENSION jchimin(6), jchimax(6), chistart(6), thetaval(ntheta), chivals(nchi)
!-----------------------------------------------------------------------
!  Determine the smallest and largest chi index for each of the
!  major divisions.
!-----------------------------------------------------------------------
nchidiv = 3
jchimin(1) = 1
jchimax(1) = 1 + 2*nlc(1)*2**ndiv
chistart(1) = zero
DO k = 2, 2*nchidiv
   IF(k <= nchidiv)kk=k
   IF(k > nchidiv)kk=k-nchidiv
   kkm1 = k-1
   IF(k > nchidiv+1)kkm1=k-nchidiv-1
   jchimin(k) = jchimax(k-1)
   jchimax(k) = jchimin(k) + 2*nlc(kk)*2**ndiv
   chistart(k) = chistart(k-1) + dchi(kkm1)*(jchimax(k-1) - jchimin(k-1) )
ENDDO
!-----------------------------------------------------------------------
!  Calculate theta and chi values.
!-----------------------------------------------------------------------
DO itheta = 1, ntheta
   theta = (itheta-1)*dthe
   thetaval(itheta) = theta
ENDDO
DO jchi = 1, nchi
   DO k = 1, 2*nchidiv
      IF(jchi >= jchimin(k) .AND. jchi <= jchimax(k)) kthdiv = k
   ENDDO
   kk = kthdiv
   IF(kthdiv > nchidiv) kk = kthdiv-nchidiv
   chi = (jchi-jchimin(kthdiv))*dchi(kk)+chistart(kthdiv)
   chivals(jchi) = chi
ENDDO
RETURN
ENDSUBROUTINE angles
