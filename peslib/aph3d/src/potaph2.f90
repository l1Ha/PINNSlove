SUBROUTINE potaph2 (rho, theta, xchi, pot, nkeyfr, kfrcnt)
USE FileUnits_Module
USE Numbers_Module
USE Numbers_Module
USE Masses_Module
USE MassFactor2_Module
USE Numeric_Kinds_Module
USE Narran_Module
!            P U R P O S E   O F    S U B R O U T I N E
! This routine calculates interparticle distances in Bohr and then
! calls surface to determine the interaction potential in Hartree's.
!            I N P U T    A R G U M E N T S
!      rho    hyperradius in Bohr.          (0 <  rho  < infinity)
!      theta  hyperangle theta in radians.  (0 < theta < pi/2)
!      xchi   hyperangle chi in radians.    (0 <  chi  < 2pi)
!            O U T P U T    A R G U M E N T S
!      pot    potential energy V(rho,theta,chi) in Hartree atomic units.
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
IMPLICIT NONE
INTEGER :: k, nkeyfr, kfrcnt, kn, keyframe, lnear (3)
REAL(Kind=WP_Kind) :: pot, rho, theta, xchi
LOGICAL :: motgraph
REAL(Kind=WP_Kind) :: mu, eta, rhox, &
 rhoz, qbig, qlit, sbig (narran), fac (narran), den (narran), &
 cthetabig (narran), sthetabig (narran), sb (narran), cb (narran), &
 carrier
real :: x (3), y (3), z (3), H (3), P (3), B (3), xscale (3), &
 yscale (3), zscale (3), tension (3), continuity (3), bias (3)

!            I N T R I N S I C    F U N C T I O N S

INTRINSIC sqrt, sin, cos
!-----------------------------------------------------------------------
!  Calculate the interparticle distances and store them in the array r.
!-----------------------------------------------------------------------
motgraph = .false.
!  Don't forget to change .bound!
IF(motgraph)THEN
   IF(kfrcnt==1)THEN
      DO kn = 1, 3
      WRITE(LW_Motion_Unit + kn, * ) 'LWMO'
      WRITE(LW_Motion_Unit + kn, * ) 1
      WRITE(LW_Motion_Unit + kn, * ) 9
      WRITE(LW_Motion_Unit + kn, * ) nkeyfr
      ENDDO
   ENDIF
   mtot = mass(1) + mass(2) + mass(3)

   mu = sqrt (mass(1) * mass(2) * mass(3) / mtot)
   d (1) = sqrt ( (mass(1) / mu) * (1 - (mass(1) / mtot) ) )
   d (2) = sqrt ( (mass(2) / mu) * (1 - (mass(2) / mtot) ) )

   d (3) = sqrt ( (mass(3) / mu) * (1 - (mass(3) / mtot) ) )
   eta = fortyfive - (theta / Two)
   rhox = rho * cos (eta)

   rhoz = rho * sin (eta)
   qbig = rhox
   qlit = rhoz

ENDIF
chi3 (1) = xchi
chi3 (2) = chi3 (1) - chij (2, 1)

chi3 (3) = chi3 (1) - chij (3, 1)
DO k = 1, narran
   s (k) = rho * sqrt (Half - Half * sin (theta) * cos (Two * chi3 (k) ) )
   sbig (k) = rho * sqrt (Half + Half * sin (theta) * cos (Two  * chi3 (k) ) )
   r (k) = d (k) * s (k)
   IF(motgraph)THEN
      fac (k) = sin (theta) * cos (Two * chi3 (k) )
      den (k) = sqrt (1 - (fac (k) * fac (k) ) )
      cthetabig (k) = (sin (theta) * sin (Two * chi3 (k) ) )   / den (k)

      sthetabig (k) = sqrt (1 - (cthetabig (k) * cthetabig (k) ) )
      sb (k) = (sthetabig (k) * s (k) * sin (Two* chi3 (k) ) ) / qbig
      carrier = ( (s (k) * sin (Two * chi3 (k) ) ) * cthetabig ( k) )
      cb (k) = (carrier + (sbig (k) * cos (Two * chi3 (k) ) ) )  / qbig

   ENDIF

ENDDO
IF(motgraph)THEN
   DO kn = 1, 3
   x (kn) = sbig (kn) * cb (kn) / d (kn)
   y (kn) = 0

   z (kn) = sbig (kn) * sb (kn) / d (kn)
   H (kn) = (360 / nkeyfr) * (kfrcnt)
   P (kn) = 0

   B (kn) = 0
   xscale (kn) = 1
   yscale (kn) = 1

   zscale (kn) = 1

   keyframe = (kfrcnt - 1) * 30
   lnear (kn) = 0
   tension (kn) = 0
   continuity (kn) = 0

   bias (kn) = 0
   WRITE(LW_Motion_Unit + kn, * ) x (kn), y (kn), z (kn), H (kn), P (kn), B (kn)
   WRITE(LW_Motion_Unit + kn, * ) xscale (kn), yscale (kn), zscale (kn)
   WRITE(LW_Motion_Unit + kn, * ) keyframe
   WRITE(LW_Motion_Unit + kn, * ) lnear (kn), tension (kn), continuity (kn),  bias (kn)
   ENDDO

ENDIF
!-----------------------------------------------------------------------
! This CALL returns the interaction potential in Hartree atomic units.
!-----------------------------------------------------------------------

CALL surface (pot, r)
RETURN
ENDSUBROUTINE potaph2
