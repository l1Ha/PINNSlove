SUBROUTINE TSI_Integrate(u, z, w, n, rstart, rend, nstpl, zed, zp,      &
                         usys2, w2, zedz, EDeriv, zeta, omega, new,   &
                         g0, g1, g2, dmu, psi, gamma,                  &
                         rhoval, xksq, energy, w0, w1, w2b,            &
                         diag, cor, asym, ndim, ipiv, ui, work)
!-----------------------------------------------------------------------
! TwoState_Inner: inner-region coupled-channel propagator (first version).
!
! Same argument list as Manolopous/logderph1 so that it can be called at
! the same place in propph1.  Instead of reconstructing the interaction
! matrix from dumped data, the production routines
!
!       potmatn(...)  +  pot_shift(...)
!
! are called at every internal substep of the sector, exactly as the
! production integrators do, so the mathematical problem is identical to
! the chain's; only the discretization differs (one-step log-derivative
! propagator, validated to machine precision in the 2x2 prototypes).
!
!   G(rmid) = (I - h^2 W/12)^-1 (I + 5 h^2 W/12)
!   Y(rmid+h) = ( G(rmid+h) - ( G(rmid) + h Y(rmid) )^-1 ) / h
!
! with W the full (unpacked) interaction matrix in the code's own units
! (potmatn output with pot_shift applied).  Y is the log-derivative
! matrix in the convention of the chain (the quantity handed between
! sectors, cf. propph1's R_Mat).
!
! Switch: APH3D_TSI_USIGN = +1 (default) uses W = +U, -1 uses W = -U,
! used to settle the sign convention of potmatn empirically.
!-----------------------------------------------------------------------
USE Numeric_Kinds_Module
USE FileUnits_Module
IMPLICIT NONE
INTEGER, PARAMETER :: lwork2 = 300
INTEGER n, nstpl, ndim, ipiv(n)
REAL(Kind=WP_Kind) usys2, rstart, rend
LOGICAL EDeriv
REAL(Kind=WP_Kind) u(*), z(n,n), w(n,n), zed(n,n), zp(n,n), w2(n,n)
REAL(Kind=WP_Kind) zedz(n,n)
REAL(Kind=WP_Kind) g0(n,n), g1(n,n), g2(n,n), dmu(n,n), psi(n), gamma(n)
REAL(Kind=WP_Kind) rhoval(*), xksq(n), energy(n), w0(n,n), w1(n,n), w2b(n,n)
REAL(Kind=WP_Kind) zeta(n), omega(n), new(n), diag(n,n), cor(n,n), asym(n,n)
REAL(Kind=WP_Kind) ui(n,n), work(n)

REAL(Kind=WP_Kind) :: Wfull(lwork2,lwork2)
REAL(Kind=WP_Kind) :: Gcur(lwork2,lwork2), Gnew(lwork2,lwork2), Tmat(lwork2,lwork2)
REAL(Kind=WP_Kind) :: Ycur(lwork2,lwork2), Ainv(lwork2,lwork2), wl(8*lwork2)
REAL(Kind=WP_Kind) :: rmid, hstep, kmax, wmax
REAL(Kind=WP_Kind) :: usign
INTEGER :: i, j, ijp, istep, nstep, iinfo, ipiv2(lwork2)
LOGICAL :: tsi_first, tsi_dump
EXTERNAL potmatn, pot_shift, DSYSV
CHARACTER(LEN=64) :: tsi_env
INTEGER :: tsi_len, tsi_stat
DATA tsi_first/.TRUE./
SAVE tsi_first, usign, tsi_dump

IF(tsi_first)THEN
   tsi_first = .FALSE.
   usign = 1.d0
   tsi_dump = .FALSE.
   CALL get_environment_variable('APH3D_TSI_USIGN', tsi_env, tsi_len, tsi_stat)
   IF(tsi_stat==0 .AND. tsi_len>0) READ(tsi_env(1:tsi_len),*) usign
   CALL get_environment_variable('APH3D_TSI_DUMP', tsi_env, tsi_len, tsi_stat)
   IF(tsi_stat==0 .AND. tsi_len>0)THEN
      IF(tsi_env(1:1)=='1' .OR. tsi_env(1:1)=='T' .OR. tsi_env(1:1)=='t') tsi_dump=.TRUE.
   ENDIF
   WRITE(Out_Unit,'(A,F6.1,A,L2)') 'TSI_Integrate: usign=', usign, '  dump=', tsi_dump
ENDIF
IF(n > lwork2) STOP 'TSI_Integrate: n > lwork2'

!-----------------------------------------------------------------------
! internal substep grid from the largest |W| at the sector start
!-----------------------------------------------------------------------
rmid = rstart
CALL potmatn (n, rmid, u, g0, g1, g2, dmu, psi, gamma, rhoval, xksq, energy, &
              w0, w1, w2b, zp, diag, cor, asym, ndim)
CALL pot_shift (n, rmid, u)

! W in the code's units; potmatn's array is packed (upper triangle, row-wise)
wmax = 0.d0
DO i = 1, n
   ijp = i*(i+1)/2
   wmax = MAX(wmax, ABS(usign*u(ijp)))
ENDDO
kmax = SQRT(MAX(wmax, 1.d-30))
hstep = 0.30d0/kmax
hstep = MIN(hstep, 5.d-3)
nstep = MAX(1, NINT((rend-rstart)/hstep))
hstep = (rend-rstart)/REAL(nstep, WP_Kind)

Ycur(1:n,1:n) = z(1:n,1:n)
CALL TSI_G(n, Wloc_of(rmid), hstep, Gcur)

DO istep = 1, nstep
   rmid = rstart + hstep*REAL(istep, WP_Kind)
   CALL potmatn (n, rmid, u, g0, g1, g2, dmu, psi, gamma, rhoval, xksq, energy, &
                 w0, w1, w2b, zp, diag, cor, asym, ndim)
   CALL pot_shift (n, rmid, u)
   CALL TSI_G(n, Wloc_of(rmid), hstep, Gnew)
   ! Y <- ( Gnew - (Gcur + h Y)^-1 ) / h
   Tmat = 0.d0
   Tmat(1:n,1:n) = Gcur(1:n,1:n) + hstep*Ycur(1:n,1:n)
   CALL TSI_Inv(n, Tmat, Ainv, iinfo)
   IF(iinfo /= 0) THEN
      WRITE(Out_Unit,*) 'TSI_Integrate: singular recursion at rmid=', rmid
      EXIT
   ENDIF
   Ycur(1:n,1:n) = (Gnew(1:n,1:n) - Ainv(1:n,1:n))/hstep
   Gcur(1:n,1:n) = Gnew(1:n,1:n)
ENDDO

z(1:n,1:n) = Ycur(1:n,1:n)

IF(tsi_dump)THEN
   OPEN(Unit=98, File='tsi_dump.txt', Status='UNKNOWN', Position='APPEND')
   WRITE(98,'(A,2ES22.14,A,I6,A,4ES15.7)') '# TSI sector ', rstart, rend, ' n=', n, &
        '  Ydiag(1:4)=', (z(i,i), i=1,MIN(4,n))
   CLOSE(98)
ENDIF

CONTAINS

FUNCTION Wloc_of(rnow) RESULT(WW)
! full (unpacked) interaction matrix at rnow from the packed potmatn output
! (this is the very array the production integrators use; potmatn/pot_shift
! were just called at rnow)
IMPLICIT NONE
REAL(Kind=WP_Kind), INTENT(IN) :: rnow
REAL(Kind=WP_Kind) :: WW(n,n)
INTEGER :: ii, jj, kk
DO ii = 1, n
   DO jj = 1, n
      IF(ii <= jj) THEN
         kk = ii + jj*(jj-1)/2
      ELSE
         kk = jj + ii*(ii-1)/2
      ENDIF
      WW(ii,jj) = usign*u(kk)
   ENDDO
ENDDO
END FUNCTION Wloc_of

SUBROUTINE TSI_G(nn, Wm, hh, Gout)
IMPLICIT NONE
EXTERNAL DSYSV
INTEGER, INTENT(IN) :: nn
REAL(Kind=WP_Kind), INTENT(IN) :: Wm(nn,nn), hh
REAL(Kind=WP_Kind), INTENT(OUT) :: Gout(nn,nn)
REAL(Kind=WP_Kind) :: A2(lwork2,lwork2), wll(8*lwork2)
INTEGER :: ipv(lwork2), inf, ii
A2(1:nn,1:nn) = -1.d0*hh*hh/12.d0*Wm(1:nn,1:nn)
Gout(1:nn,1:nn) = 5.d0*hh*hh/12.d0*Wm(1:nn,1:nn)
DO ii = 1, nn
   A2(ii,ii) = A2(ii,ii) + 1.d0
   Gout(ii,ii) = Gout(ii,ii) + 1.d0
ENDDO
CALL DSYSV('U', nn, nn, A2, lwork2, ipv, Gout, lwork2, wll, SIZE(wll), inf)
IF(inf/=0) STOP 'TSI_G: DSYSV failed'
END SUBROUTINE TSI_G

SUBROUTINE TSI_Inv(nn, Am, Aout, inf)
IMPLICIT NONE
EXTERNAL DSYSV
INTEGER, INTENT(IN) :: nn
REAL(Kind=WP_Kind), INTENT(IN) :: Am(lwork2,lwork2)
REAL(Kind=WP_Kind), INTENT(OUT) :: Aout(lwork2,lwork2)
INTEGER, INTENT(OUT) :: inf
REAL(Kind=WP_Kind) :: workm(lwork2,lwork2), wll(8*lwork2)
INTEGER :: ipv(lwork2), ii
workm(1:nn,1:nn) = Am(1:nn,1:nn)
Aout = 0.d0
DO ii = 1, nn; Aout(ii,ii) = 1.d0; ENDDO
CALL DSYSV('U', nn, nn, workm, lwork2, ipv, Aout, lwork2, wll, SIZE(wll), inf)
END SUBROUTINE TSI_Inv

END SUBROUTINE TSI_Integrate
