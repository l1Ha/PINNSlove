SUBROUTINE bitint(ntheta, nchi, idt, lnod, symmetry, jmax, ndimen)
USE FEM_Module
!-------------------------------------------------------------------
!     Routine to convert the list of turned on nodes from
!     a bit list to an INTEGER list; at the same time, it extends the list
!     to the upper half by symmetry
!-------------------------------------------------------------------

!     jmax is the maximum number of nodes in the chi direction in the lower
!     half of the range
!     idt is the INTEGER list of the possible nodes: a zero means the node
!     is being used, a -3 means the node is not being used
!     lnod is the bit list used to actually make the mesh
!-------------------------------------------------------------------


!     I N P U T    A R G U M E N T S

!     ntheta is the maximum number of nodes in the theta direction
!     nchi is the maximum number of nodes in the chi direction
!     idt      Identifies each node as free, fixed, constrained, unique, nonexistant
!     lnod
!     symmetry=false, 0<=chi<=pi ; symmetry=true, 0<=chi<=pi/2)
!     jmax
!     ndimen

IMPLICIT NONE
LOGICAL symmetry
INTEGER ntheta, nchi, jmax, ndimen, nodnum, i, j, indx, ibit, jmask
INTEGER idt(ntheta,2*jmax-1), lnod(ndimen)
INTEGER, PARAMETER:: numsto=30

nchi = 2*jmax-1
DO i=1,ntheta
   DO j=1,jmax
      ! nodnum is the position in the array of the current possible nodes
      ! indx give which word of lnod the possible node is in
      ! ibit gives which bit in the word specifies this possible node
      nodnum = ntheta*(j-1)+i
      indx = (nodnum-1)/numsto + 1
      ibit = nodnum - (indx-1)*numsto - 1
      jmask = 2**ibit
      IF(ior(lnod(indx),jmask)==lnod(indx))THEN
         idt(i,j) = Free    !     the node is turned on
      ELSE
         idt(i,j) = nonxst  !     the node is not turned on
      ENDIF
      !  extend the list by symmetry to the upper half of the allowed space
      IF(.NOT.symmetry)THEN
         idt(i,j+jmax-1) = idt(i,j)
      ELSE
         idt(i,nchi-j+1) = idt(i,j)
      ENDIF
   ENDDO
   IF(.NOT.symmetry.AND.idt(i,1)==0)idt(i,jmax) = 0
ENDDO

RETURN
END
