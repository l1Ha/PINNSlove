SUBROUTINE propph1(nt, R_Mat, w, vecnew, tstore, rbegin,          &
                    rmax, nsteps, RE_Mat, zp, usys2, w2,          &
                    zedz, eredev, EDeriv, iprt, zeta,             &
                    omega, new, g0, g1, g2, w0, w1, w2b,          &
                    psi, dmu, energy, gamma,                      &
                    xksq, ncall,smx_time,mcall, sge_time,         &
                    method, diag, cor, asym, ndim, ipiv, ui, work)
USE FileUnits_MODULE
USE Regins_Module   ! PORT (2026-10-02): for endaph
USE dip_MODULE
USE junk_MODULE
USE Boundary_Module
USE BasisRegions_Module
!
!

!             P U R P O S E    O F    S U B R O U T I N E

! Integrates second-order coupled differential equations across a
! region.
! written by G. A. Parker
!
!   Modified by M. Braunstein to treat photodissociation - Spring, 1991
!
!   Stefano's modifications of 9/11/97 added 9/12/97:
!      The variables ipiv, ui, work are now passed in as arguments,
!      which are now dimensioned by the passed-in argument nt.
!    9/16/97  TRP  Added ipiv, ui, work to argument list in CALL
!                     to Manolopous().

!             I N P U T   A R G U M E N T S

!   nt        Number of coupled second-order differential equations.
!   R_Mat      LOG derivative matrix at the begining of the region.
!   w         Temporary storage array.
!   vecnew    Temporary storage array.
!   tstore    Temporary storage array.
!   nsteps    The number of steps in the log-derivative region.
!   rbegin    Beginning distance for the log-derivative region.
!   rmax      Ending distance for the log-derivative region.
!
!   zeta      Gamma vector (inhomogeneous operator) at beginning of
!             region
!

!             O U T P U T     A R G U M E N T S

!   R_Mat    Wigner R-matrix at the end of the region.
!
!   zeta    Gamma vector (inhomogeneous operator) at end of region
!

! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

IMPLICIT NONE
CHARACTER(LEN=6) iprt
CHARACTER(LEN=10) method
LOGICAL EDeriv, ifirst
INTEGER nt, isect, nrho, nsteps, isector, i, ncall, mcall, ndim, info
INTEGER ipiv(nt)
REAL(Kind=WP_Kind) rbegin, rmax, rright, rleft, rhov1, xktot, rmax_save, rmax_orig
REAL(Kind=WP_Kind) rhov2, usys2, eredev, rhov1g, rhov2g, smx_time
REAL(Kind=WP_Kind) t_smx1, t_smx2, sge_time, t_sge1, t_sge2
REAL(Kind=WP_Kind) work(nt)
!-----------------------------------------------------------------------
! arrays dimensioned as matrices (nt*nt).
!-----------------------------------------------------------------------
REAL(Kind=WP_Kind) w(nt,nt), R_Mat(nt,nt), vecnew(nt), tstore(nt,nt)
REAL(Kind=WP_Kind) RE_Mat(nt,nt), zp(nt,nt), w2(nt,nt), zedz(nt,nt)
REAL(Kind=WP_Kind) g0(nt,nt), g1(nt,nt), g2(nt,nt), dmu(nt,nt), psi(nt)
REAL(Kind=WP_Kind) gamma(nt), rhoval(21), xksq(nt), energy(nt)
REAL(Kind=WP_Kind) ui(nt,nt), w0(nt,nt), w1(nt,nt), w2b(nt,nt)
REAL(Kind=WP_Kind) zeta(nt), omega(nt), new(nt), diag(nt,nt), cor(nt,nt)
REAL(Kind=WP_Kind) asym(nt,nt)
SAVE
COMMON /utotal/ xktot
DATA ifirst/.true./
! PORT (2026-10-02): gated dump of the log-derivative (R) matrix in the
! APH sector basis (APH3D_RDUMP / APH3D_RDUMP_LO / APH3D_RDUMP_HI).
LOGICAL rdump_on, rdump_first
REAL(Kind=WP_Kind) rdump_lo, rdump_hi
CHARACTER(LEN=64) rdump_env
INTEGER rdump_len, rdump_stat, rdump_i, rdump_j
DATA rdump_first/.true./
! PORT (2026-10-02): 2x2 inner propagator check (APH3D_TSI_CHECK=1):
! re-integrate every sector with the independent one-step propagator
! (TSI_Integrate) from the same starting matrix and compare.
LOGICAL tsi_check, tsi_init
REAL(Kind=WP_Kind) R_save(300,300), R_chain(300,300)
DATA tsi_init/.true./
! PORT (2026-10-02): 2x2 data path (block-diagonal 2N channels)
LOGICAL ps2x2, ps2_init
REAL(Kind=WP_Kind) T2(300,300)
INTEGER i2, j2, n1_2
DATA ps2_init/.true./

EXTERNAL popt, mxoutd, trnsfm, sector, overlap
!iprt='little'

IF(rdump_first)THEN
   rdump_first = .false.
   rdump_on = .false.; rdump_lo = -1.d10; rdump_hi = 1.d10
   CALL get_environment_variable('APH3D_RDUMP', rdump_env, rdump_len, rdump_stat)
   IF(rdump_stat==0 .AND. rdump_len>0)THEN
      IF(rdump_env(1:1)=='1' .OR. rdump_env(1:1)=='T' .OR. rdump_env(1:1)=='t') rdump_on=.true.
   ENDIF
   CALL get_environment_variable('APH3D_RDUMP_LO', rdump_env, rdump_len, rdump_stat)
   IF(rdump_stat==0 .AND. rdump_len>0) READ(rdump_env(1:rdump_len),*) rdump_lo
   CALL get_environment_variable('APH3D_RDUMP_HI', rdump_env, rdump_len, rdump_stat)
   IF(rdump_stat==0 .AND. rdump_len>0) READ(rdump_env(1:rdump_len),*) rdump_hi
ENDIF
IF(tsi_init)THEN
   tsi_init = .false.
   tsi_check = .false.
   CALL get_environment_variable('APH3D_TSI_CHECK', rdump_env, rdump_len, rdump_stat)
   IF(rdump_stat==0 .AND. rdump_len>0)THEN
      IF(rdump_env(1:1)=='1' .OR. rdump_env(1:1)=='T' .OR. rdump_env(1:1)=='t') tsi_check=.true.
   ENDIF
   IF(tsi_check) WRITE(Out_Unit,*) 'APH3D_TSI_CHECK active (TwoState_Inner comparison)'
   IF(tsi_check) WRITE(*,*) 'DBG propph1 init, tsi_check=', tsi_check
   IF(ps2_init)THEN
      ps2_init = .false.
      CALL Get2x2Switch(ps2x2)
      IF(ps2x2) WRITE(Out_Unit,*) 'propph1: 2x2 data path active (nt=', nt, ')'
   ENDIF
ENDIF
t_smx1 = 0.d0
t_smx2 = 0.d0
!-----------------------------------------------------------------------
! Write out intgration parameters.
!-----------------------------------------------------------------------
WRITE(Out_Unit,'(A,I5,A,ES23.15,A,ES23.15)') 'nsteps=',nsteps, ' rbegin=', rbegin, ' rmax=',rmax

!-----------------------------------------------------------------------
! calculate lower and upper bounds for energy-labeled
! channels which are to be propagated
!-----------------------------------------------------------------------
WRITE(Out_Unit,*) 'nt=',nt, '   xktot=', xktot
! PORT (2026-10-02): the energy-independent storage (pmat_unit) covers the
! APH region only.  With a Delves region requested (endaph < finish) the per
! energy path would jump its last sector from the APH end to `finish` and run
! out of storage records (APHReadn EOF detected).  End the per-energy
! propagation at endaph; the Delves region is propagated by the classic path
! (aph3d1 -> delprop -> hybridph) which writes BinOut/R_Matrix_Delves.bin.
rmax_orig = rmax
rmax_save = MIN(rmax, endaph)
rleft = rbegin
rright = rmax
rhov1 = rbegin
isect = 1
READ(Ovr_Unit)ISect,Rhov1,Rhov2
nrho = NSectors
Rhov1=Basis_Dist(Isect)
Rhov2=Basis_Dist(Isect+1)
Rleft=Sector_Boundaries(Isect)
Rright=Sector_Boundaries(Isect+1)
Rmax=Sector_Boundaries(Nsectors+1)
Nrho=Nsectors
IF(lphoto)THEN
   READ(fph_unit)rhov1g,rhov2g
ENDIF

DO isector = 1, NSectors
   IF(isector/=1)newsect=.true.
   IF(isector/=1)THEN
      READ(Ovr_Unit)ISect,Rhov1,rhov2
      Rhov1=Basis_Dist(Isector)
      Rhov2=Basis_Dist(Isector+1)
      Rleft=Sector_Boundaries(Isector)
      Rright=Sector_Boundaries(Isector+1)
      Rmax=Sector_Boundaries(Nsectors+1)
      IF(isector==NSectors)THEN
         RRight=rmax_save
         Rmax=rmax_save
      ENDIF
      Nrho=Nsectors
      IF(lphoto)THEN
         READ(fph_unit)rhov1g,rhov2g
      ENDIF
      nrho = NSectors
   ENDIF
   !-----------------------------------------------------------------------
   ! integrate from rleft to rright using the log-derivative integrator.
   !-----------------------------------------------------------------------
   IF(ifirst.AND.(iprt/='little'))THEN
       WRITE(Out_Unit,'(A,es23.15,A,es23.15,A,I5)')'rleft=',rleft,' rright=',rright,' nsteps=',nsteps
   ENDIF
   IF(tsi_check) WRITE(*,*) 'DBG sector start isector=', isector, ' tsi_check=', tsi_check
   IF(tsi_check) R_save(1:nt,1:nt) = R_Mat(1:nt,1:nt)
   IF(method=='new_logder')THEN
      !----------------------------------------------------------------------
      ! This method does not use matrix inversions but uses matrix solutions.
      !----------------------------------------------------------------------
      IF(isector==1)THEN
         WRITE(Out_Unit,'(A)')'Using new_logder'
      ENDIF
      CALL logderph1(tstore, R_Mat, w, nt, rleft, rright, nsteps, RE_Mat,  &
                     zp, usys2, w2, zedz, EDeriv, zeta, omega, new,        &
                     g0, g1, g2, dmu, psi, gamma,                          &
                     rhoval, xksq, energy, w0, w1, w2b,                    &
                     diag, cor, asym, ndim, ipiv, ui, work)
   ELSEIF(method=='man_method')THEN
      !----------------------------------------------------------------------
      ! This is the Manolopous Procedure
      !----------------------------------------------------------------------
      IF(isector==1)THEN
         WRITE(Out_Unit,*)'Using Manolopous Method'
         WRITE(Msg_Unit,*)'Using Manolopous Method'
      ENDIF
      CALL Manolopous(tstore, R_Mat, w, nt, rleft, rright, nsteps, RE_Mat,  &
                     zp, usys2, w2, zedz, EDeriv, zeta, omega, new,         &
                     g0, g1, g2, dmu, psi, gamma,                           &
                     rhoval, xksq, energy, w0, w1, w2b,                     &
                     diag, cor, asym, ndim, ipiv, ui, work)
   ELSE
      WRITE(Msg_Unit,*)'Error variable method not set correctly'
      WRITE(Msg_Unit,*)method
   ENDIF

   !-----------------------------------------------------------------------
   ! PORT (2026-10-02): 2x2 inner propagator check (TwoState_Inner).
   ! Re-integrate this sector from the same starting matrix with the
   ! independent one-step log-derivative propagator, which obtains W by
   ! calling potmatn+pot_shift itself (identical mathematical problem,
   ! different discretization), and report the difference.
   !-----------------------------------------------------------------------
   IF(tsi_check)THEN
      R_chain(1:nt,1:nt) = R_Mat(1:nt,1:nt)
      R_Mat(1:nt,1:nt) = R_save(1:nt,1:nt)
      CALL TSI_Integrate(tstore, R_Mat, w, nt, rleft, rright, nsteps, RE_Mat,  &
                         zp, usys2, w2, zedz, EDeriv, zeta, omega, new,        &
                         g0, g1, g2, dmu, psi, gamma,                          &
                         rhoval, xksq, energy, w0, w1, w2b,                    &
                         diag, cor, asym, ndim, ipiv, ui, work)
      WRITE(*,'(A,I4,A,2ES20.10,A,ES14.6,A,ES14.6)') 'TSI sector', isector, &
           ' rho=', rleft, rright, &
           '  max|dY| =', MAXVAL(ABS(R_Mat(1:nt,1:nt)-R_chain(1:nt,1:nt))), &
           '  scale=', MAXVAL(ABS(R_chain(1:nt,1:nt)))
      OPEN(Unit=95, File='tsi_check_cmp.txt', Status='UNKNOWN', Position='APPEND')
      WRITE(95,'(A,I4,A,2ES20.10,A,ES14.6,A,ES14.6,A,4ES15.7)') '# sector', isector, &
           ' rho=', rleft, rright, &
           '  max|dY| =', MAXVAL(ABS(R_Mat(1:nt,1:nt)-R_chain(1:nt,1:nt))), &
           '  scale=', MAXVAL(ABS(R_chain(1:nt,1:nt))), &
           '  chainYdiag(1:4)=', (R_chain(rdump_i,rdump_i), rdump_i=1,MIN(4,nt))
      CLOSE(95)
      R_Mat(1:nt,1:nt) = R_chain(1:nt,1:nt)
   ENDIF

   IF(iprt=='medium')THEN
      WRITE(Out_Unit,*) 'logder_matrix at r=', rright
      CALL MxOut(R_Mat, nt, nt)
      IF(EDeriv)THEN
         WRITE(Out_Unit,*) 'Derivative of logder_matrix', ' with respect to E at r=', rright
         CALL MxOut(RE_Mat, nt, nt)
      ENDIF
   ELSEIF(iprt=='full  ')THEN
      WRITE(Out_Unit,*) 'logder_matrix at r=', rright
      CALL MxOutD(R_Mat, nt, nt, 0)
      IF(EDeriv)THEN
         WRITE(Out_Unit,*) 'Derivative of logder_matrix', ' with respect to E at r=', rright
         CALL MxOutD(RE_Mat, nt, nt, 0)
      ENDIF
   ENDIF
   !-----------------------------------------------------------------------
   ! If sector adiabatic basis functions are used we need to obtain
   ! the overlap matrix and transform to the new basis.
   !-----------------------------------------------------------------------
   IF(iprt=='medium')THEN
      WRITE(Out_Unit,'(1x,A,ES15.7,A,ES15.7,A,ES15.7)')'RLeft=',RLeft,'  RRight=',RRight, ' Ovlp(1,1)=',TStore(1,1)
   ENDIF
   IF(isector /= nrho)THEN
      tstore=0.d0
      IF(ps2x2)THEN
         ! lift the single-state overlap to the block-diagonal 2N form
         n1_2 = nt/2
         CALL overlap2(rhov1, rhov2, tstore, vecnew, n1_2, ifirst)
         T2 = 0.d0
         DO i2 = 1, n1_2
            DO j2 = 1, n1_2
               T2(i2, j2)             = tstore(i2, j2)
               T2(n1_2+i2, n1_2+j2)   = tstore(i2, j2)
            ENDDO
         ENDDO
         w(1:nt,1:nt) = MATMUL(R_Mat(1:nt,1:nt), T2(1:nt,1:nt))
         R_Mat(1:nt,1:nt) = MATMUL(TRANSPOSE(T2(1:nt,1:nt)), w(1:nt,1:nt))
      ELSE
         CALL overlap2(rhov1, rhov2, tstore, vecnew, nt, ifirst)
         IF(iprt=='medium')THEN
            WRITE(Out_Unit,*)'Overlap Matrix'
            CALL MxOut(tstore, nt, nt)
         ENDIF
         w=MATMUL(R_Mat,TStore)
         R_Mat=MATMUL(TRANSPOSE(TStore),w)
      ENDIF
      IF(rdump_on .AND. rright>=rdump_lo .AND. rright<=rdump_hi)THEN
         OPEN(Unit=95, File='rdump_ph1.txt', Status='UNKNOWN', Position='APPEND')
         WRITE(95,'(A,ES22.14,A,ES22.14,A,ES22.14,A,I6)') '# R at rleft=',rleft, &
              ' rhov1=',rhov1,' rright=',rright,' nt=',nt
         DO rdump_i = 1, nt
            WRITE(95,'(500ES22.14)') (R_Mat(rdump_i,rdump_j), rdump_j = 1, nt)
         ENDDO
         WRITE(95,'(A,ES22.14,A,I6)') '# T (SF overlap) at rhov1=',rhov1,' nt=',nt
         DO rdump_i = 1, nt
            WRITE(95,'(500ES22.14)') (TStore(rdump_i,rdump_j), rdump_j = 1, nt)
         ENDDO
         CLOSE(95)
      ENDIF
      IF(EDeriv)THEN
         w=MATMUL(RE_Mat,TStore)
         RE_Mat=MATMUL(TRANSPOSE(TStore),w)
      ENDIF
      !------------------------------------------------------------------
      ! transfrom gamma - half transfrom of vector
      !------------------------------------------------------------------
      IF(lphoto)THEN
         IF(iprt=='medium')THEN
            WRITE(Out_Unit,*)'gamma vector'
            WRITE(Out_Unit,*)(zeta(i),i=1,nt)          ! WRITE intermediate gamma to file
         ENDIF  
         CALL dgemv('T',nt,nt,1.d0,tstore,nt,zeta,1,0.d0,w,1)
         DO i=1,nt
            zeta(i)=w(i,1)
         ENDDO
      ENDIF
      IF(iprt=='medium')THEN
         WRITE(Out_Unit,*) 'logder_matrix at r=', rright
         WRITE(Out_Unit,*)'transformed with an overlap',  ' matrix calculated',   &
                          ' between functions at rho = ',  rhov1,' and ',rhov2
         CALL MxOut(R_Mat, nt, nt)
         IF(EDeriv)THEN
            WRITE(Out_Unit,*) 'Derivative of logder_matrix',  ' with respect to E at r=', rright
            WRITE(Out_Unit,*)'transformed with an overlap',  ' matrix calculated',  &
                             'between functions at rho = ',  rhov1,' and ',rhov2
            CALL MxOut(RE_Mat, nt, nt)
         ENDIF
      ELSEIF(iprt=='full  ')THEN
         WRITE(Out_Unit,*) 'logder_matrix at r=', rright
         WRITE(Out_Unit,*)'transformed with an overlap',  ' matrix calculated',  &
                          'between functions at rho = ',  rhov1,' and ',rhov2
         CALL MxOutD(R_Mat, nt, nt, 0)
         IF(EDeriv)THEN
            WRITE(Out_Unit,*) 'Derivative of logder_matrix',  ' with respect to E at r=', rright
            WRITE(Out_Unit,*)'transformed with an overlap',  ' matrix calculated', &
                            'between functions at rho = ',  rhov1,' and ',rhov2
            CALL MxOutD(RE_Mat, nt, nt, 0)
         ENDIF
      ENDIF
   ENDIF
ENDDO

!-----------------------------------------------------------------------
! calculate R = Y**(-1)
!-----------------------------------------------------------------------
CALL cputime(t_smx1)
CALL unit (nt, ui)
CALL dsysv ('l', nt, nt, R_Mat, nt, ipiv, ui, nt, work, nt, info)
R_Mat=ui

CALL cputime(t_smx2)
ncall = ncall + 1
smx_time = smx_time + t_smx2-t_smx1
!-----------------------------------------------------------------------
! calculate RE = - R * YE * R
!-----------------------------------------------------------------------
IF(EDeriv)THEN
   RE_Mat=MATMUL(R_Mat,RE_Mat)
   RE_Mat=-MATMUL(RE_Mat,R_Mat)
ENDIF

rmax = rmax_orig
ifirst = .false.
RETURN
END
