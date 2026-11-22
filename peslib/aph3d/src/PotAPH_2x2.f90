SUBROUTINE PotAPH_2x2(rho, theta, xchi, p11, p22, p12)
USE Numeric_Kinds_Module
USE Numbers_Module
USE Parms_Module
USE FileUnits_Module
USE Masses_Module
USE MassFactor2_Module

!            P U R P O S E   O F    S U B R O U T I N E
! Two-electronic-state (2x2) analogue of PotAPH.  Calculates the
! interparticle distances in Bohr and returns the three independent
! elements of the symmetric 2x2 diabatic potential matrix (Hartree):
!
!            | p11  p12 |
!            | p12  p22 |
!
! p12 is odd under exchange of the two identical atoms of the diatomic
! (chi -> -chi), so an angular basis used with these matrix elements
! must be able to represent odd functions of chi.

!            I N P U T    A R G U M E N T S
!      rho    hyperradius in Bohr.          (0 <  rho  < infinity)
!      theta  hyperangle theta in radians.  (0 < theta < pi/2)
!      xchi   hyperangle chi in radians.    (0 <  chi  < 2pi)

!            O U T P U T    A R G U M E N T S
!      p11, p22, p12   2x2 potential matrix elements (Hartree).

IMPLICIT NONE

!            I N T E G E R S
INTEGER k

!            R E A L S
REAL(Kind=WP_Kind) p11, p22, p12, va1, va2, rho, theta, xchi, sin, cos

!            L O C A L    S C R A T C H
!  Local arrays instead of the like-named module arrays chi3/s/r of
!  MassFactor2_Module: those are shared state, which would make this routine
!  unusable from OpenMP threads.  The computed values are identical and no
!  caller reads them from the module.
REAL(Kind=WP_Kind) chi3l(narran), sl(narran), rl(narran)

!            I N T R I N S I C    F U N C T I O N S
INTRINSIC sqrt, sin, cos

!            E X T E R N A L S
EXTERNAL surface_krb2_2x2

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
! 2x2 diabatic potential matrix (lambda=1) in Hartree atomic units.
!-----------------------------------------------------------------------
CALL surface_krb2_2x2(rl, p11, p22, p12, va1, va2, 1.0d0)

RETURN
ENDSUBROUTINE PotAPH_2x2
