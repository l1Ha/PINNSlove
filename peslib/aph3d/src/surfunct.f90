SUBROUTINE surfunct (naph, fbar, h, idelm, phi, numnp, nquad)
USE Numeric_Kinds_Module
!
!   $RCSfile: surfunct.f,v $ $Revision: 1.15 $
!   $Date: 89/11/16 11:03:32 $
!   $State: Stable $
!

!              P U R P O S E   O F   S U B R O U T I N E

!      To evaluate the surface functions at the quadrature points
!      within a particular element.

!              I N P U T    A R G U M E N T S

!      naph         Number of surface functions.
!      h            Interpolation coefficients.
!      idelm        Nodes numbers within the element.
!      phi          Surface functions at the nodal points.
!      numnp        Number of nodal points.
!      nquad        Number of Gauss_Legendre quadrature points.

!              O U T P U T    A R G U M E N T S

!      fbar         Surface functions evaluated at the Gauss_Legendre
!                   quadrature points.

IMPLICIT NONE
INTEGER i, j, idelm, imap, k, n, naph, numnp, nquad
REAL(Kind=WP_Kind) fbar, h, phi, zero
DIMENSION  fbar(naph, nquad, nquad), h(9, nquad, nquad),  phi(numnp, naph), imap(9), idelm(9)
EXTERNAL vsets

PARAMETER (zero=0.d0)
!              D A T A    S T A T E M E N T S

DATA imap / 2, 6, 3, 5, 9, 7, 1, 8, 4 /

!-----------------------------------------------------------------------
!  Zero the surface function array.
!-----------------------------------------------------------------------
CALL vsets (naph*nquad*nquad,  fbar, 1, zero)
!-----------------------------------------------------------------------
! Use the interpolation coefficients to evaluate the surface functions
! at the quadrature points.
!-----------------------------------------------------------------------
DO 85 j = 1, nquad
   DO 83 i = 1, nquad
      DO 80 k = 1, 9
         DO 79 n = 1,naph
            fbar(n, i, j) = fbar(n, i, j)+h(imap(k), i, j) *phi(idelm(k), n)
   79          CONTINUE
   80       CONTINUE
   83    CONTINUE
   85 CONTINUE
RETURN
!
END
