MODULE Das_Module
!-----------------------------------------------------------------------
!  dscale  Mass scaling factors for obtaining mass weighted coordinates.
!          S = dscale * R
!          s = r / dscale
!          where R and r the unscaled Jacobi coordinates and S and s
!          are the scaled or mass weighted Jacobi coordinates.
!          Coordinate r is the distance between the particles in the
!          diatomic fragment and R is the distance between the center
!          of mass of the diatomic fragment and the remaining particle.
!-----------------------------------------------------------------------
   USE Numeric_Kinds_Module
   USE Parms_Module
   IMPLICIT NONE
   SAVE
   REAL(Kind=WP_Kind) dscale(narran) !  dscale  Mass scaling factors, S = dscale * R,  s = r / dscale.
END MODULE Das_Module
