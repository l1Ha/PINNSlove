SUBROUTINE transf(umat, raph, rdelv, naph, ndelv)
!---------------------------------------------------------------------------------
!Modified by: Xuan LI, Department of Physics and Astronomy, University of Oklahoma
!---------------------------------------------------------------------------------

!            P U R P O S E   O F    S U B R O U T I N E

!            I N P U T    A R G U M E N T S
!  umat
!  raph
!  rdelv
!  naph
!  ndelv

!            O U T P U T    A R G U M E N T S
!This routine is called by:
!                          main_states
!
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
USE fileunits_Module
USE Numeric_Kinds_Module
IMPLICIT NONE
INTEGER naph, ndelv
REAL(Kind=WP_Kind) umat(ndelv, naph),raph(naph, naph), rdelv(ndelv, ndelv)

rdelv=Matmul(Matmul(umat,raph),Transpose(umat))

RETURN
ENDSUBROUTINE transf
