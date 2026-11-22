SUBROUTINE thetachi (th, ch, thetaval, chivals, ntheta, nchi, numnp, idelm, thetas, chis)
USE Numeric_Kinds_Module
!
!   $RCSfile: thetachi.f,v $ $Revision: 1.15 $
!   $Date: 89/11/16 11:03:34 $
!   $State: Stable $
!


!              P U R P O S E   O F   S U B R O U T I N E

!  Returns the theta and chi coordinates of the nodal points.


!                  I N P U T   A R G U M E N T S

!   th        List of theta indicies for the theta direction.
!   ch        List of chi indicies for the chi direction.
!   thetaval  List of theta values.
!   chivals   List of chi values.
!   ntheta    Number of theta angles.
!   nchi      Number of chi angles.
!   numnp     Number of nodal points.
!   idelm     List of node numbers.


!               O U T P U T   A R G U M E N T S

!   thetas    Value of theta at each node within the element.
!   chis      Value of chi at each node within the element.

IMPLICIT NONE
INTEGER k, th, ch, numnp, ntheta, nchi, idelm, node
REAL(Kind=WP_Kind) thetaval, chivals, thetas, chis
DIMENSION th(numnp), ch(numnp), thetaval(ntheta), chivals(nchi), thetas(9), chis(9), idelm(9)
!-----------------------------------------------------------------------
!  Determine theta and chi values.
!-----------------------------------------------------------------------

DO k = 1,9
   node = idelm(k)
   thetas(k) = thetaval(th(node))
   chis(k)   = chivals(ch(node))
ENDDO

RETURN
!
END
