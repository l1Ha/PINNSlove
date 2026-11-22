SUBROUTINE surface_krb2(v, rij)
!-----------------------------------------------------------------------
! PORT (2026-10-01): K + Rb2 / Rb + KRb2 potential for the APH3D port.
!
! APH3D interface:  rij(1:3) = internuclear distances (bohr) where rij(k) is
! the distance of the diatom pair of arrangement channel k, i.e. the pair
! that excludes atom k (see BoundDiatomic/PotDiatomic usage).
!
! Atom labelling assumed by the input deck:  1 = K (lightest), 2 = Rb,
! 3 = Rb.  Hence
!     channel 1 (K spectator)    : rij(1) = Rb-Rb
!     channel 2 (Rb spectator)   : rij(2) = K-Rb
!     channel 3 (Rb spectator)   : rij(3) = K-Rb
!
! KRb2 ESI PES convention (krb2_pes_eval): q(1)=K1-Rb, q(2)=K2-Rb, q(3)=Rb-Rb,
!   giving the 2x2 diabatic matrix [[Vd11,Vd12],[-Vd12,Vd22]] and the adiabatic
!   eigenvalues Va1 <= Va2 (hartree), with lambda = 1 (diabatic) ... 0 (adiabatic).
!
! This routine returns the LOWER ADIABATIC surface (the single-surface
! "adiabatic" calculations of the reference paper).  The 2x2 diabatic matrix
! is available through surface_krb2_2x2 (used by the multi-state extension).
!-----------------------------------------------------------------------
USE Numeric_Kinds_Module
USE Masses_Module
IMPLICIT NONE
REAL(Kind=WP_Kind), INTENT(OUT) :: v
REAL(Kind=WP_Kind), INTENT(IN)  :: rij(3)
REAL*8 :: q(3), vd11, vd22, vd12, Va1, Va2
LOGICAL, SAVE :: krb2_ready = .FALSE.
INTEGER :: i
REAL*8 :: wall_ev
REAL*8, PARAMETER :: ev_to_ha = 1.0d0/27.211386245988d0

IF (.NOT. krb2_ready) THEN
   CALL krb2_pes_init
   krb2_ready = .TRUE.
ENDIF

! map APH3D channel distances to the PES (K-Rb1, K-Rb2, Rb-Rb) convention
q(1) = rij(3)    ! K - Rb(atom 2)
q(2) = rij(2)    ! K - Rb(atom 3)
q(3) = rij(1)    ! Rb - Rb
! Numerical guards for the RKHS fit: the kernel overflows (=> NaN) for
! R >~ 1e25 (pothyp passes rtrans = 1e30 as "infinity") and extrapolates
! poorly below ~0.5 bohr.  Clamping at 1e6 bohr is exact to machine
! precision since all long-range terms have vanished there.
q = MIN(MAX(q, 0.5d0), 1.0d6)

! lambda = 1: full physical PES = 2x2 diabatic matrix (pairwise diatomic
! part + RKHS three-body interaction; the latter is switched off smoothly
! for R > 28 bohr by SW).  Va1 is the lower eigenvalue of that matrix, i.e.
! the adiabatic ground surface used by the single-surface chain.
! (lambda = 0 would scale out the three-body interaction and leave only the
! pairwise diatomic part - that is NOT the physical potential.)
CALL krb2_pes_eval(q, 1.0d0, vd11, vd22, vd12, Va1, Va2)

! Short-range repulsive wall: the RKHS fit extrapolates poorly (goes to ~0)
! below ~5 bohr; same compensation as in the ABC project's amod_pes_multi
! (10 eV * (5.5 - R)^4 per pair, converted to hartree).
wall_ev = 0.0d0
DO i = 1, 3
   IF (q(i) < 5.5d0) wall_ev = wall_ev + 10.0d0*((5.5d0 - q(i))**4)
ENDDO
v = Va1 + wall_ev*ev_to_ha

RETURN
END SUBROUTINE surface_krb2

SUBROUTINE surface_krb2_2x2(rij, vd11, vd22, vd12, Va1, Va2, lambda)
!-----------------------------------------------------------------------
! 2x2 diabatic matrix version of the KRb2 potential (atomic units).
!   lambda = 1 -> diabatic 2x2 representation
!   lambda = 0 -> adiabatic (Vd11 = Va1, Vd22 = Va2, Vd12 = 0)
!-----------------------------------------------------------------------
USE Numeric_Kinds_Module
USE Masses_Module
IMPLICIT NONE
REAL(Kind=WP_Kind), INTENT(IN)  :: rij(3), lambda
REAL(Kind=WP_Kind), INTENT(OUT) :: vd11, vd22, vd12, Va1, Va2
REAL*8 :: q(3)
LOGICAL, SAVE :: krb2_ready2 = .FALSE.
INTEGER :: i
REAL*8 :: wall_ev
REAL*8, PARAMETER :: ev_to_ha = 1.0d0/27.211386245988d0

IF (.NOT. krb2_ready2) THEN
   CALL krb2_pes_init
   krb2_ready2 = .TRUE.
ENDIF
q(1) = rij(3)
q(2) = rij(2)
q(3) = rij(1)
q = MIN(MAX(q, 0.5d0), 1.0d6)

IF (lambda >= 0.999d0) THEN
   CALL krb2_pes_eval(q, 1.0d0, vd11, vd22, vd12, Va1, Va2)
ELSE
   CALL krb2_pes_eval(q, 0.0d0, vd11, vd22, vd12, Va1, Va2)
ENDIF

! short-range repulsive wall (see surface_krb2)
wall_ev = 0.0d0
DO i = 1, 3
   IF (q(i) < 5.5d0) wall_ev = wall_ev + 10.0d0*((5.5d0 - q(i))**4)
ENDDO
vd11 = vd11 + wall_ev*ev_to_ha
vd22 = vd22 + wall_ev*ev_to_ha
Va1  = Va1  + wall_ev*ev_to_ha
Va2  = Va2  + wall_ev*ev_to_ha

RETURN
END SUBROUTINE surface_krb2_2x2
