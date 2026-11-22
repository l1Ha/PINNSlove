SUBROUTINE adini (id, theta, chi, nid, idi, beta, nodj, nel)
USE FileUnits_Module
USE Masses_Module
!
!   $RCSfile: adini.f,v $ $Revision: 1.15 $
!   $Date: 89/11/16 11:02:50 $
!   $State: Stable $
!

!            P U R P O S E   O F    S U B R O U T I N E
!-----------------------------------------------------------------------
! This routine reads in input DATA for the finite element method.
!-----------------------------------------------------------------------


!            I N P U T    A R G U M E N T S

!     id(numnp)   id=-nce if the is constrained and positive otherwise
!     theta
!     chi
!     nid(NDISCE) Array with the number of terms in the constraint equations
!     idi
!     beta(ndisce) Constraint Coefficients
!     nodj
!     nel


IMPLICIT NONE
LOGICAL little, medium, full
INTEGER id, nid(1), idi, nodj, nq, neq, nel, numnp, ndisce, nidm, nfreq, nmode, ithcall, ithsub
REAL(Kind=WP_Kind) beta, theta, chi, rtol
DIMENSION id(1), theta(1), chi(1), idi(1), nodj(1), beta(1)
EXTERNAL inputx, coneqn

!-----------------------------------------------------------------------
!  Determine printing options.
!-----------------------------------------------------------------------
DATA ithcall/0/, ithsub/0/
DATA little/.false./, medium/.false./, full/.false./
CALL popt ('adini   ', little, medium, full, ithcall, ithsub)

READ(Msher_Bin_Unit)  nel, numnp, ndisce, nidm, nfreq, nq, nmode, rtol
nmode=nfreq

IF(little)THEN
   WRITE(Out_unit,*) nel, numnp, ndisce, nidm, nfreq, nq, nmode, rtol
ENDIF

!-----------------------------------------------------------------------
! READ nodal point DATA.
!-----------------------------------------------------------------------
CALL inputx (id, theta, chi, numnp, neq)

!-----------------------------------------------------------------------
! READ constraint equations DATA.
!-----------------------------------------------------------------------
CALL coneqn (id, nid, idi, beta, nodj, ndisce, nidm, numnp)

RETURN
ENDSUBROUTINE Adini
