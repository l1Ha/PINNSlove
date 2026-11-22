SUBROUTINE logderph2(u, z, w, n, rstart, rend, nstpl, zed, zp,     &
                     usys2, w2, zedz, EDeriv, zeta, omega, new,    &
                     g0, g1, g2, dmu, psi, gamma,                  &
                     rhoval, xksq, energy, w0, w1, w2b,            &
                     ncall, smx_time, mcall, sge_time, diag, cor,  &
                     asym, ndim, ipiv, ui, work)

!            P U R P O S E   O F    S U B R O U T I N E
!  This routine propagates a set of coupled second-order differential
!  equations using Bob Johnson's Log Derivative method.

!            I N P U T    A R G U M E N T S
!  u       An temporary arrary for storing the interaction matrix
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
LOGICAL EDeriv
INTEGER n, nstpl, nstep, istep, nread, ncall, mcall, ndim, ii
REAL(Kind=WP_Kind) rstart, rend, r, d1, d2, d3, d4, drnow, wavel, h
REAL(Kind=WP_Kind) usys2, derinit, hstep, smx_time, t_smx1, t_smx2
REAL(Kind=WP_Kind) sge_time, t_sge1, t_sge2
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

! Calculate energy independent.
wavel=2.d0*pi/sqrt(usys2*engmax)
IF(wavel.gt.2.d0) wavel=2.d0
nstep=nstpl*(rend-rstart)/wavel+.9d0
IF(nstep.eq.0) nstep=1
drnow=(rend-rstart)/nstep
h=drnow*0.5d0
d1=h*h/3.d0
d2=2.d0*d1
d4=-d1/16.d0
d3=2.d0*d2
IF(lphoto)THEN
   READ(fph_unit)wavel,nstep,drnow,h,d1,d2,d3
ENDIF
!
! Calculate the interaction matrix.
! NOTE: The negative of the interaction matrix returned.
!
r=rstart
CALL potmatn ( n, r, u, g0, g1, g2, dmu, psi, gamma, rhoval, xksq, energy, &
               w0, w1, w2b, zp, diag, cor, asym, ndim)

CALL pot_shift (n, r, u)
!
IF(lphoto.and.r.lt.enddip)THEN
   CALL phot_READ(n, omega)
ELSEIF(lphoto)THEN
   CALL vsets(n, omega,1,Zero)
ENDIF
!
! determine z(0) +I = h*y(0) - h**2*u(0)/3 + I
!
CALL vscale(n*n,       z,   1, h )
CALL vscale(n*(n+1)/2, u,   1, d1)
CALL specadd (n, z, u)
CALL onead (n, z)
IF(EDeriv)THEN
   CALL vscale(n*n,       zed, 1, h )
   CALL addhone (n, zed, n, -d1*usys2)
ENDIF
!
! determine zeta(0)=new(0) + h*gamma(0)
!
IF(lphoto)THEN
   CALL vscale(n, zeta,1,h)
   CALL vscale(n, omega,1,d1)
   CALL vadd(n, zeta, omega)
ENDIF
!
! start integrating the differential equations.
!
DO istep = 1, nstep
   !
   ! determine the value of r.
   !
   r=rstart+h*(2*istep-1)
   !
   ! determine the potential matrix elements
   !
   CALL potmatn ( n, r, u, g0, g1, g2, dmu, psi, gamma, rhoval, xksq, energy,  &
                  w0, w1, w2b, zp, diag, cor, asym, ndim)
   CALL pot_shift (n, r, u)
   !
   IF(lphoto.and.r.le.enddip)THEN
      CALL phot_READ(n, omega)
   ELSEIF(lphoto)THEN
      CALL vsets(n, omega, 1, Zero)
   ENDIF
   !
   ! determine [1/8+h**2*u(2*n+1)/48]
   !
   CALL vscale(n*(n+1)/2, u, 1, d4)
   CALL vsets (n*n, w, 1, Zero)
   CALL specadd (n, w, u)
   CALL addhone (n, w, n, One/Eight)
   !
   ! determine [1/8+h**2*u(2*n+1)/48]**(-1)
   !
   CALL cputime(t_smx1)
   CALL unit (n, ui)
   CALL dsysv ('l', n, n, w, n, ipiv, ui, n, work, n, info)
   CALL cputime(t_smx2)
   ncall = ncall + 1
   smx_time = smx_time + t_smx2-t_smx1
   DO ii = 1,n
     CALL dcopy (n, ui(1,ii), 1, w(1,ii), 1)
   ENDDO
   !
   ! determine new(1)=(I-(h^2/6)*u(1))^-1 * (4h^2/3)*omega(1)
   !
   IF(lphoto)THEN
      CALL vscale(n, omega, 1, d1)
      CALL vscale(n, omega, 1, Half)
      CALL dgemv('N',n,n,1.d0,w,n,omega,1,0.d0,new,1)
   ENDIF
   !
   ! invert I+z
   !
   CALL cputime(t_smx1)
   CALL unit (n, ui)
   CALL dsysv ('l', n, n, z, n, ipiv, ui, n, work, n, info)
   DO ii = 1,n
      CALL dcopy (n, ui(1,ii), 1, z(1,ii), 1)
   ENDDO
   CALL cputime(t_smx2)
   ncall = ncall + 1
   smx_time = smx_time + t_smx2-t_smx1

   CALL fill_mat(z,n)
   !
   ! determine zeta(1)=new(1) + Z(0)^-1 * zeta(0)
   !
   IF(lphoto)THEN
      CALL dgemv('N',n,n,1.d0,z,n,zeta,1,Zero,omega,1)
      CALL vadd (n, new, omega)
      CALL vsets (n, zeta, 1, Zero)
      CALL vadd (n, zeta, new)
   ENDIF
!
!   determine the energy derivative of (z + i)(2*n + 1)
!
   IF(EDeriv)THEN
      CALL cputime(t_sge1)
      CALL dgemm('n', 'n', n,n,n,1d0,w,n,w,n,Zero, w2, n)
      CALL cputime(t_sge2)
      sge_time = sge_time + t_sge2-t_sge1
      CALL vscale (n*n, w2, 1, -d4*usys2)
      CALL cputime(t_sge1)
      CALL dgemm('t', 'n', n, n, n,1d0,zed,n,z,n,Zero, zedz, n)
      CALL dgemm('t', 'n', n,n,n,1d0,z,n,zedz,n,Zero, zed, n)
      sge_time = sge_time + t_sge2-t_sge1
      mcall = mcall + 3
      CALL vscale(n*n, zed, 1, -One)
      CALL vadd (n*n, zed, w2)
   ENDIF
   !
   ! determine [I+z(2*n+1)]=[1/8+h**2*u(2*n+1)/48]**(-1)
   !                    -[I+z(2*n)]**(-1)-6*I
   !
   CALL vscale(n*n, z, 1, -One)
   CALL vadd (n*n, z, w)
   CALL addhone(n, z, n, -Six)
   !***********************************************************************
   !                   Start determining (z + i)(2n* + 2)
   !                   and its energy derivative
   !***********************************************************************
   !
   !
   ! determine [I+z(2*n+1)]**(-1)
   !
   CALL cputime(t_smx1)
   CALL unit (n, ui)
   CALL dsysv ('l', n, n, z, n, ipiv, ui, n, work, n, info)
   DO ii = 1,n
      CALL dcopy (n, ui(1,ii), 1, z(1,ii), 1)
   ENDDO
   CALL cputime(t_smx2)
   ncall = ncall + 1
   smx_time = smx_time + t_smx2-t_smx1
   !
   ! determine the value of r.
   !
   r=rstart+h*(2*istep)
   !
   ! determine the potential matrix elements
   !
   CALL potmatn ( n, r, u, g0, g1, g2, dmu, psi, gamma,  rhoval, xksq, energy, &
                  w0, w1, w2b, zp, diag, cor, asym, ndim)
   CALL pot_shift (n, r, u)
   !
   IF(lphoto.and.r.le.enddip)THEN
      CALL phot_READ(n, omega)
   ELSEIF(lphoto)THEN
      CALL vsets(n, omega, 1, Zero)
   ENDIF
   !
   IF(istep.eq.nstep) d2=d1
   !
   !  determine the energy derivative of z(2*n + 2)
   !
   IF(EDeriv)THEN
      CALL cputime(t_sge1)
      CALL dgemm('t', 'n', n,n,n,1d0,zed,n,z,n,Zero,zedz, n)
      CALL dgemm('t', 'n', n,n,n,1d0,z,n,zedz,n,Zero, zed, n)
      CALL cputime(t_sge2)
      sge_time = sge_time + t_sge2-t_sge1
      mcall = mcall + 2
      CALL vscale (n*n, zed, 1, -One)
      CALL addhone(n, zed, n, -d2*usys2)
   ENDIF
   !
   ! determine zeta(2)=new(2) + Z(1)^-1 * zeta(1)
   !
   IF(lphoto)THEN
      CALL vscale (n, omega, 1, d2)
      CALL vsets (n, new, 1, Zero)
      CALL vadd (n, new, omega)  
      CALL dgemv('N',n,n,One,z,n,zeta,1,Zero,omega,1)
      CALL vadd (n, new, omega)
      CALL vsets(n, zeta, 1, Zero)
      CALL vadd (n, zeta, new)
   ENDIF
   !
   !
   ! determine [z(2*n + 2)+I]=-[I+z(2*n-1)]**(-1)-2*h**2*u(2*n)/3+2I
   !
   CALL vscale (n*n, z, 1, -One)
   CALL vscale (n*(n+1)/2, u, 1, d2)
   CALL specadd(n, z, u)
   CALL addhone(n, z, n, Two)
ENDDO
!
! determine y(2*nstep)=[z(2*n)+I]/h-I/h
!
CALL vscale (n*n, z, 1, One/h)
CALL addhone (n, z, n, -One/h)
!
! detrmine gamma_2 at end of sector
!
IF(lphoto)THEN
   CALL vscale(n, zeta, 1, One/h)
ENDIF

IF(EDeriv)THEN
    CALL vscale (n*n, zed, 1, One/h)
ENDIF
!***********************************************************************
!     now z = Y, zed = YE
!***********************************************************************
RETURN
END
