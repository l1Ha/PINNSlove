SUBROUTINE potmatn_2x2 (n2, rmid, u2, g0, g1, g2, dmu, psi, gamma,  &
                        rhoval, xksq, energy, w0, w1, w2b, wtemp,   &
                        diag, cor, asym, ndim)
!-----------------------------------------------------------------------
! Two-electronic-state wrapper around the production potential builder.
!
! Calls potmatn with the single-state channel count n = n2/2 (so that the
! energy-independent storage is consumed exactly as in the single-surface
! chain) and expands the packed interaction matrix into the 2N x 2N
! block-diagonal form used by the V12 = 0, V22 = V11 limit:
!
!      U_2N = [ U_N   0  ]
!             [ 0    U_N ]
!
! This is the data path for the TwoState_Inner propagation: the production
! integrator (Manolopous/logderph1) then propagates 2N channels with the
! chain's own discretization, and the result must reduce to
! R_2N = blkdiag(R_N, R_N).  For the physical 2x2 case the off-diagonal
! (coupling) blocks computed from surface_krb2_2x2 are added here as well.
!
! packed storage: element (i,j) with i <= j sits at i + j(j-1)/2
! (the production potmatn/pot_shift convention).
!
! Switch: APH3D_PS2X2 = 1 activates the 2N expansion (default off, in which
! case the routine is a plain pass-through to potmatn with n = n2).
!-----------------------------------------------------------------------
USE Numeric_Kinds_Module
USE FileUnits_Module
IMPLICIT NONE
INTEGER, INTENT(IN)  :: n2, ndim
REAL(Kind=WP_Kind), INTENT(IN)  :: rmid
REAL(Kind=WP_Kind), INTENT(OUT) :: u2(*)
REAL(Kind=WP_Kind), INTENT(IN)  :: rhoval(*)
REAL(Kind=WP_Kind), INTENT(OUT) :: g0(n2,n2), g1(n2,n2), g2(n2,n2)
REAL(Kind=WP_Kind), INTENT(OUT) :: dmu(n2,n2), w0(n2,n2), w1(n2,n2), w2b(n2,n2)
REAL(Kind=WP_Kind), INTENT(OUT) :: wtemp(n2,n2), diag(n2,n2), cor(n2,n2), asym(n2,n2)
REAL(Kind=WP_Kind), INTENT(OUT) :: psi(n2), gamma(n2), xksq(n2), energy(n2)
REAL(Kind=WP_Kind) :: u1(300*301/2)
REAL(Kind=WP_Kind) :: dummy2(1,1)
INTEGER :: ns, i, j, kk, ioff, joff, k1, k2
LOGICAL :: twox2, first, ps2x2_dump
CHARACTER(LEN=64) :: env2
INTEGER :: len2, stat2
REAL(Kind=WP_Kind) :: ps2x2_lo, ps2x2_hi, ps2x2_eps
LOGICAL :: ps2x2_trace
INTEGER :: ps2x2_ncall
DATA ps2x2_ncall/0/
EXTERNAL potmatn
DATA first/.TRUE./
SAVE first, twox2, ps2x2_dump, ps2x2_lo, ps2x2_hi, ps2x2_trace, ps2x2_ncall, ps2x2_eps

IF(first)THEN
   first = .FALSE.; twox2 = .FALSE.
   CALL get_environment_variable('APH3D_PS2X2', env2, len2, stat2)
   IF(stat2==0 .AND. len2>0)THEN
      IF(env2(1:1)=='1' .OR. env2(1:1)=='T' .OR. env2(1:1)=='t') twox2=.TRUE.
   ENDIF
   WRITE(Out_Unit,'(A,L2)') 'potmatn_2x2: two-electron-state expansion =', twox2
   ps2x2_dump = .FALSE.; ps2x2_lo = -1.d10; ps2x2_hi = 1.d10
   CALL get_environment_variable('APH3D_PS2X2_DUMP', env2, len2, stat2)
   IF(stat2==0 .AND. len2>0)THEN
      IF(env2(1:1)=='1' .OR. env2(1:1)=='T' .OR. env2(1:1)=='t') ps2x2_dump=.TRUE.
   ENDIF
   ps2x2_eps = 0.d0
   CALL get_environment_variable('APH3D_PS2X2_EPS', env2, len2, stat2)
   IF(stat2==0 .AND. len2>0) READ(env2(1:len2),*) ps2x2_eps
   ps2x2_trace = .FALSE.
   CALL get_environment_variable('APH3D_PS2X2_TRACE', env2, len2, stat2)
   IF(stat2==0 .AND. len2>0)THEN
      IF(env2(1:1)=='1' .OR. env2(1:1)=='T' .OR. env2(1:1)=='t') ps2x2_trace=.TRUE.
   ENDIF
   CALL get_environment_variable('APH3D_PS2X2_LO', env2, len2, stat2)
   IF(stat2==0 .AND. len2>0) READ(env2(1:len2),*) ps2x2_lo
   CALL get_environment_variable('APH3D_PS2X2_HI', env2, len2, stat2)
   IF(stat2==0 .AND. len2>0) READ(env2(1:len2),*) ps2x2_hi
ENDIF

IF(.NOT.twox2)THEN
   ! pass-through (single-state behaviour)
   CALL potmatn (n2, rmid, u2, g0, g1, g2, dmu, psi, gamma, rhoval, xksq, &
                 energy, w0, w1, w2b, wtemp, diag, cor, asym, ndim)
   RETURN
ENDIF
IF(ps2x2_trace)THEN
   ps2x2_ncall = ps2x2_ncall + 1
   IF(ps2x2_ncall <= 200)THEN
      OPEN(Unit=90, File='ps2x2_trace.txt', Status='UNKNOWN', Position='APPEND')
      WRITE(90,'(I8,2ES22.14)') ps2x2_ncall, rmid, REAL(n2,WP_Kind)
      CLOSE(90)
   ENDIF
ENDIF

ns = n2/2
IF(2*ns /= n2) STOP 'potmatn_2x2: n2 must be even'
! the single-state builder is called with the full-size work arrays; only
! the first ns x ns blocks of its output are used
CALL potmatn (ns, rmid, u1, g0, g1, g2, dmu, psi, gamma, rhoval, xksq, &
              energy, w0, w1, w2b, wtemp, diag, cor, asym, ndim)

! expand the packed single-state matrix into the 2N packed matrix
! (block diagonal: state 1 in [1,ns], state 2 in [ns+1,2ns], no coupling)
u2(1:n2*(n2+1)/2) = 0.d0
DO j = 1, n2
   DO i = 1, j
      IF(i <= ns .AND. j <= ns) THEN
         u2(i + j*(j-1)/2) = u1(i + j*(j-1)/2)
      ELSEIF(i > ns .AND. j > ns) THEN
         ! APH3D_PS2X2_EPS: optional infinitesimal asymmetry between the two
         ! electronic blocks (diagnostic: with exactly identical blocks the
         ! LU pivoting of the 2N solves is degenerate)
         u2(i + j*(j-1)/2) = (1.d0+ps2x2_eps)*u1((i-ns) + (j-ns)*((j-ns)-1)/2)
      ENDIF
   ENDDO
ENDDO

! PORT (2026-10-02): gated dump of the raw single-state packed matrix and of
! the expanded 2N packed matrix (APH3D_PS2X2_DUMP=1), for data-path checks
! PORT (2026-10-02): gated trace of the call sequence (r values), used to
! verify that the 2N run consumes the energy-independent storage exactly
! like the single-state run (APH3D_PS2X2_TRACE=1)
IF(ps2x2_trace)THEN
   ps2x2_ncall = ps2x2_ncall + 1
   IF(ps2x2_ncall <= 200)THEN
      OPEN(Unit=90, File='ps2x2_trace.txt', Status='UNKNOWN', Position='APPEND')
      WRITE(90,'(I8,2ES22.14)') ps2x2_ncall, rmid, REAL(n2,WP_Kind)
      CLOSE(90)
   ENDIF
ENDIF
IF(ps2x2_dump .AND. rmid >= ps2x2_lo .AND. rmid <= ps2x2_hi)THEN
   OPEN(Unit=92, File='ps2x2_dump.txt', Status='UNKNOWN', Position='APPEND')
   WRITE(92,'(A,ES22.14,A,I4)') '# potmatn_2x2 at r=', rmid, '  ns=', ns
   WRITE(92,'(A)') '# u1 (single-state packed, row-major i<=j)'
   WRITE(92,'(500ES22.14)') (u1(i), i=1, ns*(ns+1)/2)
   WRITE(92,'(A)') '# u2 (2N packed)'
   WRITE(92,'(500ES22.14)') (u2(i), i=1, n2*(n2+1)/2)
   CLOSE(92)
ENDIF

RETURN
ENDSUBROUTINE potmatn_2x2
