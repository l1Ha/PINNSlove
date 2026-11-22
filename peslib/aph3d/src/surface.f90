SUBROUTINE surface(v, rij)
!-----------------------------------------------------------------------
! PORT (2026-10-01, gfortran port of APH3D):
!   System dispatcher for the APH3D potential interface.  In the original
!   distribution each system shipped its own surface.f90 that both checked
!   PES_Name and called the system potential; here the check is central and
!   the system routine is selected by PES_Name.
!
!   On entry:  rij(1:3) = internuclear distances in bohr (not mass scaled).
!              rij(k) is the distance of the diatom pair of arrangement
!              channel k, i.e. the pair excluding atom k
!              (see BoundDiatomic: varying rij(k) with the others at
!               infinity gives channel k's diatomic curve).
!   On exit:   v = potential energy in hartree.
!-----------------------------------------------------------------------
USE PES_Module
USE Numeric_Kinds_Module
IMPLICIT NONE
REAL(Kind=WP_Kind), INTENT(OUT) :: v
REAL(Kind=WP_Kind), INTENT(IN)  :: rij(3)
CHARACTER(LEN=25) :: name

name = PES_Name

IF (INDEX(name,'DMBE') > 0 .OR. INDEX(name,'PKH3') > 0 .OR. INDEX(name,'Porter-Karplus') > 0) THEN
   CALL surface_dmbe(v, rij)
ELSEIF (INDEX(name,'KRb2') > 0 .OR. INDEX(name,'KRB2') > 0 .OR. INDEX(name,'KRb') > 0) THEN
   CALL surface_krb2(v, rij)
ELSE
   WRITE(0,*) 'surface: unknown PES_Name: ', name
   STOP 'surface: unknown PES_Name'
ENDIF
RETURN
END SUBROUTINE surface
