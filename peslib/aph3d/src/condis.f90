SUBROUTINE condis (nid, idi, beta, phi, nidm, ndisce, neq)
USE Numeric_Kinds_Module
!-----------------------------------------------------------------------
! This routine determines the values of the surface functions at the
! nodal points that have been constrained to be a linear combination
! of the values of the corresponding surface function at unconstrained
! nodal points.
!-----------------------------------------------------------------------

!  I N P U T    A R G U M E N T S

!  nid         number of constraints for nodal point i.
!  idi
!  beta        constraint coefficient.
!  phi         surface function at nodal point k.
!  nidm
!  ndisce      number of nodal points that are constrained.
!  neq


IMPLICIT NONE
INTEGER nid, idi, nidm, ndisce, i, neq, j
REAL(Kind=WP_Kind) beta, phi, phik, zero
DIMENSION nid(1), idi(nidm,1), beta(nidm,1), phi(1)
PARAMETER (zero= 0.d0)
!-----------------------------------------------------------------------
!   phi(k)    surface function at nodal point k.
!   ndisce    number of nodal points that are constrained.
!   beta(j)   constraint coefficient.
!   nid(i)    number of constraints for nodal point i.
!-----------------------------------------------------------------------

DO i=1, ndisce
   phik= zero
   DO j=1, nid(i)
      phik=phik+beta(j, i)*phi(idi(j, i))
   ENDDO
   phi(neq+i)=phik
ENDDO

RETURN
ENDSUBROUTINE condis
