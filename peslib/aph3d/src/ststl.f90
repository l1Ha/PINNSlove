SUBROUTINE ststl(c, theta, phi)
USE FileUnits_Module
USE rhosur_Module
USE totj_Module
USE pass_Module
!----------------------------------------------------------------------
!  This routine calculates the coefficients (c) that are used in the
!  variational equation.
!----------------------------------------------------------------------
IMPLICIT NONE
LOGICAL full
INTEGER i, j
REAL(Kind=WP_Kind) c, theta, phi, wfunc, gtheta, gphi, zero, two, four
DIMENSION c(4,4)
EXTERNAL poteff, gtotal
PARAMETER (zero=0.d0, two=2.d0, four=4.d0)

DATA full/.false./
!----------------------------------------------------------------------
!  Determine the interaction potential.
!  This CALL also determines the interatomic distances s(karran)
!  which are stored in common save.
!----------------------------------------------------------------------
CALL poteff (rhosurf, theta, phi, wfunc)
!----------------------------------------------------------------------
!  Calculating the weighting function.
!----------------------------------------------------------------------
CALL gtotal(theta, phi, gtot, gtheta, gphi)
!----------------------------------------------------------------------
! Initialize all coefficients to zero.
!----------------------------------------------------------------------
DO j = 1, 4
   DO i = 1, 4
      c(i,j) = zero
   ENDDO
ENDDO
!----------------------------------------------------------------------
! Set up the coefficients in the variational equation.
!----------------------------------------------------------------------
c(4,4) = wfunc*sin(two*theta)*gtot*gtot+four*sin(two*theta)*gtheta*gtheta+two *cos(theta)/sin(theta)*gphi*gphi
c(1,1) = four*sin(two*theta)*gtot*gtot
c(3,3) = two*cos(theta)/sin(theta)*gtot*gtot
c(1,4) = four*sin(two*theta)*gtheta*gtot
c(4,1) = c(1,4)
c(3,4) = two*cos(theta)/sin(theta)*gphi*gtot
c(4,3) = c(3,4)

IF(full)THEN
   WRITE(Out_unit,*)' rhosurf=',rhosurf,' theta=', theta,' phi=',phi
   WRITE(Out_unit,*)'gtot=',gtot,' gtheta=',gtheta,' gphi=',gphi
   WRITE(Out_unit,*)'  Coefficient matrix'
   DO i=1,4
      WRITE(Out_unit,7)(c(i,j),j=1,4)
   ENDDO
ENDIF
   7  FORMAT(1x,5e14.7)
RETURN
END
