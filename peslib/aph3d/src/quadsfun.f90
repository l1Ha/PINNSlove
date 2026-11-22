SUBROUTINE quadsfun (thetas, chis, hamil, fbar, fbarp, value, naph, naphp, h, p, xpt, wht, nquad)
USE Numeric_Kinds_Module
USE rhosur_Module
USE rhoovr_Module
USE covlp_Module
USE FileUnits_Module
USE Numbers_Module
!
!   $RCSfile: quadsfun.f,v $ $Revision: 1.15 $
!   $Date: 89/11/16 11:03:21 $
!   $State: Stable $
!

!              P U R P O S E   O F   S U B R O U T I N E


!                  I N P U T   A R G U M E N T S

!     thetas     List of theta values at the nodes.
!     chis       List of chi values at the nodes.
!     fbar       Surface functions at the quadrature points.
!     fbarp      Surface functions at the quadrature points which
!                need not be the same functions as fbar.
!     value      Value of the integrand at the quadrature points.
!     naph       Number of surface functions in array fbar.
!     naphp      Number of surface functions in array fbarp.
!     p          Interpolation coefficients for evaluating the
!                jacobian.
!     nquad      Number of Gauss_Legendre quadrature points.
!     xpt        Abscissas for the Gauss_Legendre quadrature.
!     wht        Weights for the Gauss_Legendre quadrature.


!               O U T P U T   A R G U M E N T S

!     hamil      hamiltonian matrix elements

IMPLICIT NONE
INTEGER i, j, imap(9), k, n, np, naph, naphp, nquad
REAL(Kind=WP_Kind) chi, theta, wij, det, gtheta, gchi,  gtot1, gtot2
REAL(Kind=WP_Kind) chis(9), thetas(9), h(9, nquad, nquad), fbar(naph, nquad, nquad)
REAL(Kind=WP_Kind) fbarp(naphp, nquad, nquad)
REAL(Kind=WP_Kind) hamil(naph, naphp), xj(2, 2), p(2, 9, nquad, nquad), value(nquad, nquad)
REAL(Kind=WP_Kind) xpt(nquad), wht(nquad)
INTRINSIC sin
EXTERNAL grms

!               D A T A    S T A T E M E N T S

DATA imap /2, 6, 3, 5, 9, 7, 1, 8, 4/

!-----------------------------------------------------------------------
! evaluate the hamiltonian matrix element.
!-----------------------------------------------------------------------
DO j = 1,nquad
   DO i = 1,nquad
      wij = wht(i)*wht(j)

!-----------------------------------------------------------------------
! Determine theta and chi at the quadrature point.
!-----------------------------------------------------------------------
      theta = zero
      chi = zero
      xj(1, 1) = zero
      xj(2, 1) = zero
      xj(1, 2) = zero
      xj(2, 2) = zero
      DO k = 1, 9
         theta = theta+thetas(k)*h(imap(k), i, j)
         chi = chi+chis(k)*h(imap(k), i, j)
         xj(1, 1) = xj(1, 1) + p(1, imap(k), i, j)*thetas(k)
         xj(2, 1) = xj(2, 1) + p(2, imap(k), i, j)*thetas(k)
         xj(1, 2) = xj(1, 2) + p(1, imap(k), i, j)*chis(k)
         xj(2, 2) = xj(2, 2) + p(2, imap(k), i, j)*chis(k)
      ENDDO

!-----------------------------------------------------------------------
! Evaluate the Jacobian.
!-----------------------------------------------------------------------
      IF(.NOT.lovlp)THEN
         CALL grms(theta, chi, gtot1, gtheta, gchi)
         gtot2 = gtot1
      ELSE
         rhosurf = rho1t
         CALL grms(theta, chi, gtot1, gtheta, gchi)
         rhosurf = rho2t
         CALL grms(theta, chi, gtot2, gtheta, gchi)
      ENDIF
      det = xj(1, 1)*xj(2, 2) - xj(2, 1)*xj(1, 2)
      wij = wij*det*gtot1*gtot2*sin(2*theta)*value(i,j)

!-----------------------------------------------------------------------
! Store matrix element contribution into the hamiltonian matrix.
!-----------------------------------------------------------------------
      IF(lovlp)THEN
         DO np = 1, naphp
            DO n = 1, naph
               hamil(n, np) = hamil(n, np)+fbar(n, i, j)*fbarp(np, i, j)*wij
            ENDDO
         ENDDO
      ELSE
         DO np = 1, naphp
            DO n = 1, np
               hamil(n, np) = hamil(n, np)+fbar(n, i, j) *fbarp(np, i, j)*wij
            ENDDO
         ENDDO
      ENDIF
   ENDDO
ENDDO

IF(.NOT.lovlp)THEN
   DO np = 1, naphp
      DO n = 1, np
         hamil(np, n) = hamil(n, np)
      ENDDO
   ENDDO
ENDIF

RETURN
END
