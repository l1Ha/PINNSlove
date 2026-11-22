SUBROUTINE coneqn (id, nid, idi, beta, nodj, ndisce, nidm, numnp)
USE Numeric_Kinds_Module
Use FileUnits_Module
!-----------------------------------------------------------------------
! This routine reads in the constraint equations and establishes constraint
! equation numbers.

!            I N P U T    A R G U M E N T S
!    numnp             Number of nodal points
!    ndisce            Number of nodes to be constrained.
!    nidm              Maximum number of terms in the constraint equation.
!    id(numnp)         id=-nce if the is constrained and positive otherwise
!
!            O U T P U T    A R G U M E N T S
!    id(numnp)         id=-nce if the is constrained and positive otherwise
!    nid(NDISCE)       Array with the number of terms in the constraint equations
!    idi(nidm,ndisce)  id=-nce if the is constrained and positive otherwise
!    beta(nidm,ndisce) Constraint Coefficients
!    nodj(nidm)        Node numbers for constraint equations
!-----------------------------------------------------------------------

IMPLICIT NONE
LOGICAL little, medium, full
INTEGER                         :: nce, j, ncecheck, nodn, jj, ithcall, ithsub
INTEGER, INTENT(IN)             :: ndisce, numnp,  nidm 
INTEGER, INTENT(OUT)            :: nid(ndisce), idi(nidm,ndisce), nodj(nidm)
INTEGER, INTENT(INOUT)          :: id(numnp)
REAL(Kind=WP_Kind), INTENT(OUT) :: beta(nidm,ndisce)
EXTERNAL popt
DATA ithcall/0/, ithsub/0/
DATA little/.false./, medium/.false./, full/.false./
CALL popt ('coneqn  ', little, medium, full, ithcall, ithsub)
!-----------------------------------------------------------------------
! Loop over the constrained equations.
!-----------------------------------------------------------------------
 DO  nce=1,ndisce
   !-----------------------------------------------------------------------
   ! READ the following:
   !   ncecheck  constraint equation number.
   !   nodn      code number to be constrained.
   !   nid(nce)  number of terms in the constraint equation.
   !-----------------------------------------------------------------------
   READ(Msher_Bin_Unit) ncecheck, nodn, nid(nce)
   IF(full) WRITE(Out_unit,*) 'coneqn: ',ncecheck, nodn, nid(nce)
   !-----------------------------------------------------------------------
   ! Make sure that nodn is less than or equal to nunmp.
   !-----------------------------------------------------------------------
   IF(nodn>numnp)THEN
      WRITE(Out_unit,*)'error in coneqn: nodn>nump ', nodn, numnp
      STOP  'coneqn'
   ENDIF
   !-----------------------------------------------------------------------
   ! Make sure that nid(nce) is less than or equal to nidm.
   !-----------------------------------------------------------------------
   IF(nid(nce)>nidm)THEN
      WRITE(Out_unit,*)'error in coneqn: nid(nce)>nidm ', nid(nce), nidm
      STOP 'coneqn'
   ENDIF
   !-----------------------------------------------------------------------
   ! Make sure that this is a constrained node.
   !-----------------------------------------------------------------------
   IF(id(nodn)/=-2)THEN
      WRITE(Out_unit,*)'error in coneqn: id(nodn)/=-2, :', id(nodn)
      STOP 'coneqn'
   ENDIF
   !-----------------------------------------------------------------------
   ! Set the id(nodn) equal to -nce which labels it as a constrained node.
   !-----------------------------------------------------------------------
   id(nodn)=-nce
   !-----------------------------------------------------------------------
   ! READ in the following:
   !   nodj        node numbers for constraint equations.
   !   beta        constraint coefficients.
   !-----------------------------------------------------------------------
   READ(Msher_Bin_Unit) (nodj(j),beta(j,nce),j=1,nid(nce))
   IF(full) WRITE(Out_unit,*) 'coneqn: ', (nodj(j),  beta(j,nce),j=1,nid(nce))
   !-----------------------------------------------------------------------
   !     establish equation number of independent degrees of freedom
   !-----------------------------------------------------------------------
   DO j=1,nid(nce)
      jj=nodj(j)
      !-----------------------------------------------------------------------
      ! Make sure that nodn is less than or equal to numnp.
      !-----------------------------------------------------------------------
      IF(nodj(j)>numnp)THEN
         WRITE(Out_unit,*)'error in coneqn: nodj(j)>numnp ', nodj(j), numnp
         STOP 'coneqn'
      ENDIF
      !-----------------------------------------------------------------------
      ! Make sure that the free nodes are indeed free.
      !-----------------------------------------------------------------------
      IF(id(jj)<=0)THEN
         WRITE(Out_unit,*)'nce,jj = ', nce,jj
         STOP 'coneqn'
      ENDIF
      idi(j,nce)=id(jj)
   ENDDO
ENDDO

RETURN
END
