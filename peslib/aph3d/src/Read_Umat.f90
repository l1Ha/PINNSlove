SUBROUTINE Read_Umat(naph, ndelve, umat, chanl, elect, nvib, jrot, lorb, energy, umat_unit)
!---------------------------------------------------------------------------------
!Modified by: Xuan LI, Department of Physics and Astronomy, University of Oklahoma
!This routine is called by:
!                          main_umatph
!---------------------------------------------------------------------------------

USE Numeric_Kinds_Module
IMPLICIT NONE
INTEGER naph, ndelve, chanl(ndelve), elect(ndelve), nvib(ndelve), jrot(ndelve), lorb(ndelve), umat_unit
REAL(Kind=WP_Kind)energy(ndelve), umat(ndelve, naph)

READ(umat_unit)energy
READ(umat_unit)chanl
READ(umat_unit)elect
READ(umat_unit)nvib
READ(umat_unit)jrot
READ(umat_unit)lorb
READ(umat_unit)umat

RETURN
ENDSUBROUTINE Read_Umat
