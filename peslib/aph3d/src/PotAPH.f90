SUBROUTINE PotAPH(rho, theta, xchi, pot)
USE Numeric_Kinds_Module
USE Numbers_Module
USE Parms_Module
USE FileUnits_Module
USE Masses_Module
USE MassFactor2_Module

!            P U R P O S E   O F    S U B R O U T I N E
! This routine calculates interparticle distances in Bohr and THEN
! calls surface to determine the interaction potential in Hartree's.

!            I N P U T    A R G U M E N T S
!      rho    hyperradius in Bohr.          (0 <  rho  < infinity)
!      theta  hyperangle theta in radians.  (0 < theta < pi/2)
!      xchi   hyperangle chi in radians.    (0 <  chi  < 2pi)

!            O U T P U T    A R G U M E N T S
!      pot    potential energy V(rho,theta,chi) in Hartree atomic units.

IMPLICIT NONE

!            I N T E G E R S
INTEGER k

!            R E A L S
REAL(Kind=WP_Kind) pot, rho, theta, xchi, sin, cos

!            L O C A L    S C R A T C H
!  Local arrays instead of the like-named module arrays chi3/s/r of
!  MassFactor2_Module: those are shared state, which would make this routine
!  unusable from OpenMP threads.  The computed values are identical and no
!  caller reads them from the module.
REAL(Kind=WP_Kind) chi3l(narran), sl(narran), rl(narran)

!            I N T R I N S I C    F U N C T I O N S
INTRINSIC sqrt, sin, cos

!            E X T E R N A L S
EXTERNAL surface

!-----------------------------------------------------------------------
!  Calculate the interparticle distances and store them in the array r.
!-----------------------------------------------------------------------
chi3l(1)=xchi
chi3l(2)=chi3l(1)-chij(2,1)
chi3l(3)=chi3l(1)-chij(3,1)
DO k=1,narran
   sl(k)=rho*sqrt(Half-Half*sin(theta)*cos(Two*chi3l(k)))
   rl(k)=d(k)*sl(k)
ENDDO

!-----------------------------------------------------------------------
! This CALL returns the interaction potential in Hartree atomic units.
!-----------------------------------------------------------------------
CALL surface (pot, rl)

RETURN
ENDSUBROUTINE PotAPH
