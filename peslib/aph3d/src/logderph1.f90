SUBROUTINE logderph1(u, z, w, n, rstart, rend, nstpl, zed, zp,     &
                     usys2, w2, zedz, EDeriv, zeta, omega, new,    &
                     g0, g1, g2, dmu, psi, gamma,                  &
                     rhoval, xksq, energy, w0, w1, w2b,            &
                     diag, cor, asym, ndim, ipiv, ui, work)

!            P U R P O S E   O F    S U B R O U T I N E
!  This routine propagates a set of coupled second-order differential
!  equations using Bob Johnson's Log Derivative method.

!            I N P U T    A R G U M E N T S
!  u       An temporary array for storing the interaction matrix
!  z       On entering contains the Log Derivative matrix at rstart.
!  n       Number of coupled channels.
!  rstart  Starting distance.
!  rend    Ending distance.
!  nstpl   Number of steps per asymptotic wavelength.
!
!  zeta    gamma vector at rstart
!  omega   temporary vector used in gamma propagation
!  new     temporary vector used in gamma propagation
!

!            O U T P U T    A R G U M E N T S
! z        On return z contains the Wiger R-matrix at rend.
!
! zeta     gamma vector at rend
!

! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!trp>>
!  Modified on 9/12/97:  Changed ui, work, and ipiv to arguments of
!     dimension n,n.  Changed calls to dsysv appropriately.  Eliminated
!     references to lw.  Eliminated some no-longer-needed diagnostic
!     output.
!  Modified on 9/15/97:  Eliminated some further no-longer-needed
!     diagnostic output.
!    9/16/97  TRP  Modified calls to addhone() in accordance with
!                  changes in that routine.
!<<trp
USE Numbers_Module
USE dip_MODULE
USE engpro_MODULE
USE FileUnits_MODULE
IMPLICIT NONE
! PORT (2026-10-02): 2x2 (two electronic state) data path switch
LOGICAL ps2x2, ps2x2_first
DATA ps2x2_first/.true./
SAVE ps2x2_first, ps2x2
EXTERNAL potmatn_2x2
LOGICAL EDeriv
INTEGER n, nstpl, nstep, istep, nread, ncall, mcall, ndim, ii
REAL(Kind=WP_Kind) rstart, rend, r, d1, d2, d3, d4, drnow, wavel, h
REAL(Kind=WP_Kind) usys2, derinit, hstep
REAL(Kind=WP_Kind) u(n,n), z(n,n), w(n,n), zed(n,n), zp(n,n), w2(n,n)
REAL(Kind=WP_Kind) zedz(n,n)
REAL(Kind=WP_Kind) g0(n,n), g1(n,n), g2(n,n), dmu(n,n), psi(n), gamma(n)
REAL(Kind=WP_Kind) rhoval(*), xksq(n), energy(n), w0(n,n), w1(n,n), w2b(n,n)
REAL(Kind=WP_Kind) zeta(n), omega(n), new(n), diag(n,n), cor(n,n), asym(n,n)
REAL(Kind=WP_Kind) ui(n,n), work(n)
INTEGER ipiv(n), info
SAVE
DATA derinit / 1.d0 /
EXTERNAL pot_shift, smxinv
EXTERNAL phot_read

! Calculate energy independent parameters
wavel=2.d0*pi/sqrt(usys2*engmax)
IF(wavel>2.d0) wavel=2.d0
nstep=nstpl*(rend-rstart)/wavel+.9d0
IF(nstep==0) nstep=1
drnow=(rend-rstart)/nstep
h=drnow*0.5d0
d1=h*h/3.d0
! PORT (2026-10-02): 2x2 (two electronic state) data path switch
IF(ps2x2_first)THEN
   ps2x2_first = .false.; ps2x2 = .false.
   CALL Get2x2Switch(ps2x2)
   WRITE(Out_Unit,'(A,L2)') 'logderph1: two-electron-state data path =', ps2x2
ENDIF
d2=2.d0*d1
d4=-d1/16.d0
d3=2.d0*d2
IF(lphoto)THEN
   READ(fph_unit)wavel,nstep,drnow,h,d1,d2,d3
ENDIF
!
! Calculate the interaction matrix.
! NOTE: The negative of the interaction matrix returned.
! NOTE: The interaction Matrix is stored in symmetric packed form.
r=rstart
! PORT (2026-10-02): two-electron-state (2x2) mode feeds the propagator
! with the block-diagonal 2N interaction matrix (see potmatn_2x2)
IF(ps2x2)THEN
   CALL potmatn_2x2 (n, r, u, g0, g1, g2, dmu, psi, gamma, rhoval, xksq, energy, &
                     w0, w1, w2b, zp, diag, cor, asym, ndim)
ELSE
   CALL potmatn ( n, r, u, g0, g1, g2, dmu, psi, gamma, rhoval, xksq, energy, &
                  w0, w1, w2b, zp, diag, cor, asym, ndim)
ENDIF

CALL pot_shift (n, r, u)
CALL Zdump1(n, z, r)
!
IF(lphoto.AND.r<enddip)THEN
   CALL phot_READ(n, omega)
ELSEIF(lphoto)THEN
   omega=0.d0
ENDIF
!
! determine z(0) +I = h*y(0) - h**2*u(0)/3 + I
!
Z=z*h
u=u*d1
CALL specadd (n, z, u)
CALL onead (n, z)
IF(EDeriv)THEN
   zed=zed*h
   CALL addhone (n, zed, n, -d1*usys2)
ENDIF
!
! determine zeta(0)=new(0) + h*gamma(0)
!
IF(lphoto)THEN
   zeta=zeta*h
   omega=omega*d1
   zeta=zeta+omega
ENDIF
!
! start integrating the differential equations.
!
DO istep = 1, nstep
   r=rstart+h*(2*istep-1)   ! determine the value of r.
   !
   ! determine the potential matrix elements
   !
   IF(ps2x2)THEN
      CALL potmatn_2x2 (n, r, u, g0, g1, g2, dmu, psi, gamma, rhoval, xksq, energy,  &
                        w0, w1, w2b, zp, diag, cor, asym, ndim)
   ELSE
      CALL potmatn ( n, r, u, g0, g1, g2, dmu, psi, gamma, rhoval, xksq, energy,  &
                     w0, w1, w2b, zp, diag, cor, asym, ndim)
   ENDIF
   CALL pot_shift (n, r, u)
   !
   IF(lphoto.AND.r<=enddip)THEN
      CALL phot_READ(n, omega)
   ELSEIF(lphoto)THEN
      omega=0.d0
   ENDIF
   !
   ! determine [1/8+h**2*u(2*n+1)/48]
   !
   u=u*d4
   w=0.d0
   CALL specadd (n, w, u)
   CALL addhone (n, w, n, One/Eight)
   !
   ! determine [1/8+h**2*u(2*n+1)/48]**(-1)
   !
   CALL unit (n, ui)
   CALL dsysv ('l', n, n, w, n, ipiv, ui, n, work, n, info)
   IF(info/=0) WRITE(*,'(A,I6,A,I6)') 'DBG logderph1 dsysv(info/=0): n=', n, ' info=', info
   w=ui
   !
   ! determine new(1)=(I-(h^2/6)*u(1))^-1 * (4h^2/3)*omega(1)
   !
   IF(lphoto)THEN
      omega=omega*d1
      omega=omega*0.5d0
      CALL dgemv('N',n,n,1.d0,w,n,omega,1,0.d0,new,1)
   ENDIF
   !
   ! invert I+z
   !
   CALL unit (n, ui)
   CALL dsysv ('l', n, n, z, n, ipiv, ui, n, work, n, info)
   IF(info/=0) WRITE(*,'(A,I6,A,I6)') 'DBG logderph1 dsysvZ(info/=0): n=', n, ' info=', info
   z=ui
   !
   ! determine zeta(1)=new(1) + Z(0)^-1 * zeta(0)
   !
   IF(lphoto)THEN
      CALL dgemv('N',n,n,1.d0,z,n,zeta,1,0.d0,omega,1)
      new=new+omega
      zeta=new
   ENDIF
!
!   determine the energy derivative of (z + i)(2*n + 1)
!
   IF(EDeriv)THEN
      CALL dgemm('n', 'n', n,n,n,1d0,w,n,w,n,0d0, w2, n)
      w2=-d4*usys2*w2
      CALL dgemm('n', 'n', n, n, n,1d0, zed,n,z,n,0d0, zedz, n) !!
      CALL dgemm('n', 'n', n, n, n,1d0,z,n,zedz,n,0d0,  zed, n) !!
      zed=w2-zed
   ENDIF
   !
   ! determine [I+z(2*n+1)]=[1/8+h**2*u(2*n+1)/48]**(-1)
   !                    -[I+z(2*n)]**(-1)-6*I
   !
   z=w-z
   CALL addhone(n, z, n, -six)
   !***********************************************************************
   !                   Start determining (z + i)(2n* + 2)
   !                   and its energy derivative
   !***********************************************************************
   !
   !
   ! determine [I+z(2*n+1)]**(-1)
   !
   CALL unit (n, ui)
   CALL dsysv ('l', n, n, z, n, ipiv, ui, n, work, n, info)
   IF(info/=0) WRITE(*,'(A,I6,A,I6)') 'DBG logderph1 dsysvZ(info/=0): n=', n, ' info=', info
   z=ui
   !
   ! determine the value of r.
   !
   r=rstart+h*(2*istep)
   !
   ! determine the potential matrix elements
   !
   IF(ps2x2)THEN
      CALL potmatn_2x2 ( n, r, u, g0, g1, g2, dmu, psi, gamma,  rhoval, xksq, energy, &
                        w0, w1, w2b, zp, diag, cor, asym, ndim)
   ELSE
      CALL potmatn ( n, r, u, g0, g1, g2, dmu, psi, gamma,  rhoval, xksq, energy, &
                     w0, w1, w2b, zp, diag, cor, asym, ndim)
   ENDIF
   CALL pot_shift (n, r, u)
   !
   IF(lphoto.AND.r<=enddip)THEN
      CALL phot_READ(n, omega)
   ELSEIF(lphoto)THEN
      omega=0.d0
   ENDIF
   !
   IF(istep==nstep) d2=d1
   !
   !  determine the energy derivative of z(2*n + 2)
   !
   IF(EDeriv)THEN
      CALL dgemm('n', 'n', n,n,n,1d0, zed,n,z,n,0d0,zedz, n)
      CALL dgemm('n', 'n', n,n,n,1d0,z,n,zedz,n,0d0, zed, n)
      zed=-1.d0*zed
      CALL addhone(n, zed, n, -d2*usys2)
   ENDIF
   !
   ! determine zeta(2)=new(2) + Z(1)^-1 * zeta(1)
   !
   IF(lphoto)THEN
      omega=omega*d2
      new=omega 
      CALL dgemv('N',n,n,1.d0,z,n,zeta,1,0.d0,omega,1)
      new=new+Omega
      zeta=new
   ENDIF
   !
   !
   ! determine [z(2*n + 2)+I]=-[I+z(2*n-1)]**(-1)-2*h**2*u(2*n)/3+2I
   !
   z=-1.d0*z
   u=u*d2
   CALL specadd(n, z, u)
   CALL addhone(n, z, n, Two)
ENDDO
!
! determine y(2*nstep)=[z(2*n)+I]/h-I/h
!
z=z/h
CALL addhone (n, z, n, -One/h)
!
! detrmine gamma_2 at end of sector
!
IF(lphoto)THEN
   zeta=zeta/h
ENDIF

IF(EDeriv)THEN
   zed=zed/h
ENDIF
!***********************************************************************
!     now z = Y, zed = YE
!***********************************************************************
RETURN
CONTAINS
SUBROUTINE Zdump1(nn, zz, rr)
! PORT (2026-10-02): gated dump of the integrator's state z
USE FileUnits_MODULE
IMPLICIT NONE
INTEGER, INTENT(IN) :: nn
REAL(Kind=WP_Kind), INTENT(IN) :: zz(nn,nn), rr
LOGICAL :: zd_on, zd_first
REAL(Kind=WP_Kind) :: zd_lo, zd_hi
CHARACTER(LEN=64) :: zd_env
INTEGER :: zd_len, zd_stat, zi
DATA zd_first/.true./
SAVE zd_on, zd_lo, zd_hi
IF(zd_first)THEN
   zd_first = .false.
   zd_on = .false.; zd_lo = -1.d10; zd_hi = 1.d10
   CALL get_environment_variable('APH3D_ZDUMP', zd_env, zd_len, zd_stat)
   IF(zd_stat==0 .AND. zd_len>0)THEN
      IF(zd_env(1:1)=='1' .OR. zd_env(1:1)=='T' .OR. zd_env(1:1)=='t') zd_on=.true.
   ENDIF
   CALL get_environment_variable('APH3D_ZDUMP_LO', zd_env, zd_len, zd_stat)
   IF(zd_stat==0 .AND. zd_len>0) READ(zd_env(1:zd_len),*) zd_lo
   CALL get_environment_variable('APH3D_ZDUMP_HI', zd_env, zd_len, zd_stat)
   IF(zd_stat==0 .AND. zd_len>0) READ(zd_env(1:zd_len),*) zd_hi
ENDIF
IF(.NOT.zd_on) RETURN
IF(rr < zd_lo .OR. rr > zd_hi) RETURN
OPEN(Unit=97, File='zdump.txt', Status='UNKNOWN', Position='APPEND')
WRITE(97,'(A,ES22.14,A,I6,A,4ES16.8)') '# z(r) at r=', rr, '  n=', nn, &
     '  diag(1:4)=', (zz(zi,zi), zi = 1, MIN(4,nn))
CLOSE(97)
END SUBROUTINE Zdump1
END
