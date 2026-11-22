SUBROUTINE boundold (nrho, rho_val, ntheta, tht_val, nchi, chi_val,    &
 drho, dtheta, dchi, psi, mval, dafx, hermit, ndaf, nchanl, nd_daf,    &
 toe, eig, v_pot, diag_opt, nsize, t_rho, t_theta, t_chi, t_chir,      &
 t_chi1, t_chir1, system, group, a, b, c, d, rhom2, jtot, usys2,       &
 direct, tblock, iblock, hamil)
USE Numeric_Kinds_Module
USE Numbers_Module
USE FileUnits_Module
USE Convrsns_Module
!-----------------------------------------------------------------------
! This routine was written by G. A. Parker
! If you find an error or have an improvement please send a messge to
! Parker@ou.edu
!-----------------------------------------------------------------------
IMPLICIT NONE  
CHARACTER(LEN=3) :: system, group, repname
CHARACTER(LEN=15) :: Identifier  
LOGICAL :: direct, debug, fulldiag  
INTEGER :: mval, jval, ndaf, kmax, nchanl, nchan, i, nd_daf,  info
INTEGER :: nbound, diag_opt, mfound, ifail, nrho, ntheta, nchi, nsize
INTEGER :: nchir, n1 (8), nirrep, irrep, nsymop, jtot, parity, par (8)
INTEGER :: itheta, irho, iblock (nchi, nchi), npot, nsize_new, lammin
INTEGER :: lammax, nbegin, fac
REAL(Kind=WP_Kind) :: dafx (0:nd_daf, 0:ndaf), abstol, prefac, eiga1,  &
 hermit (0:mval + ndaf + 1), tblock (nchi, nchi), drho, dtheta,        &
 dchi, psi (nsize, nsize), toe (0:nsize-1), eig (nsize), v_pot (nsize),&
 hamil (nsize, nsize), rho_val (nrho), tht_val (ntheta),               &
 chi_val (nchi), t_rho (nrho, nrho), t_theta (ntheta, ntheta),         &
 t_chi (nchi, nchi), t_chir (nchi, nchi), a (ntheta), b (ntheta),      &
 c (ntheta), d (ntheta), rhom2 (nrho), usys2, t_chi1 (nchi, nchi),     &
 t_chir1 (nchi, nchi), cpu_begin, cpu_end, eig_im (1000)

REAL(Kind=WP_Kind), ALLOCATABLE :: psi_new (:, :)  
WRITE(Out_Unit, * )
WRITE(Out_Unit, * ) 'Begin(BoundOld)'
debug = .False.
!-----------------------------------------------------------------------
! Calculate A, B, C, and D
!-----------------------------------------------------------------------
DO itheta = 1, ntheta
   IF(Itheta==1)THEN
     a(1)=1.d0
     b(1)=1.d0
     c(1)=1.d0
     d(1)=1.d0
   ELSE
      a (itheta) = 1.d0 / (1.d0 + sin (tht_val (itheta) ) )
      b (itheta) = 0.5d0 / sin (tht_val (itheta) ) **2
      c (itheta) = 1.d0 / (1.d0 - sin (tht_val (itheta) ) )
      d (itheta) = cos (tht_val (itheta) ) / sin (tht_val (itheta) ) **2
   ENDIF
   IF (debug)THEN
      WRITE(Out_Unit, '(1x,A15,i5,5es14.7)') 'itheta,a,b,c,d=', itheta, a(itheta) , b(itheta) , c(itheta) , d(itheta)
   ENDIF
ENDDO

b (1) = 10000.d0
!-----------------------------------------------------------------------
! Calculate 1/rho**2
!-----------------------------------------------------------------------
DO irho = 1, nrho
   IF(rho_val (irho)==0.d0)THEN
      rhom2 (irho)=Ten**(RANGE(One)-10)
   ELSE
      rhom2 (irho) = One / rho_val (irho) **2 
   ENDIF
   IF (debug)THEN
      WRITE(Out_Unit, * ) 'irho,rhom2=', irho, rhom2 (irho)
   ENDIF
ENDDO
!----------------------------------------------------------------------
! Calculate rho DAF
!----------------------------------------------------------------------
abstol = 1.d-10
IF (debug) WRITE(Out_Unit, * ) 'T_rho'
IF (nrho>1)THEN
   CALL cputime (cpu_begin)
   Identifier="_BoundOld_"
   CALL get_daf (nrho, rho_val, drho, mval, dafx, hermit, ndaf,  nd_daf, toe, usys2, kmax, Identifier)
   CALL cputime (cpu_end)
   WRITE(Out_Unit, * ) 'get_daf: time=', cpu_end-cpu_begin
ENDIF
!----------------------------------------------------------------------
! Calculate rho Kinetic energy term.
!----------------------------------------------------------------------
CALL cputime (cpu_begin)
CALL get_kinetic (nrho,   toe, t_rho,   rho_val, kmax, 'rho     ', system, group, tblock, iblock, usys2, dafx, nd_daf, ndaf)
CALL cputime (cpu_end)
WRITE(Out_Unit, * ) 'get_kinetic: time=', cpu_end-cpu_begin

!----------------------------------------------------------------------
! Calculate theta DAF
!----------------------------------------------------------------------
IF (debug) WRITE(Out_Unit, * ) 'T_theta'
IF (ntheta>1)THEN
   CALL cputime (cpu_begin)
   Identifier="_theta"
   CALL get_daf (ntheta, tht_val, dtheta, mval, dafx, hermit, ndaf, nd_daf, toe, usys2, kmax, Identifier)
   CALL cputime (cpu_end)
   WRITE(Out_Unit, * ) 'get_daf: time=', cpu_end-cpu_begin
ENDIF

!---------------------------------------------------------------------
! Calculate theta Kinetic energy term
!---------------------------------------------------------------------
CALL cputime (cpu_begin)
CALL get_kinetic (ntheta, toe, t_theta, tht_val, kmax, 'theta   ', system, group, tblock, iblock, usys2, dafx, nd_daf, ndaf)
CALL cputime (cpu_end)
WRITE(Out_Unit,*)cpu_end,cpu_begin,ntheta, "pos4"
WRITE(Out_Unit,*)"get_kinetic: time=", cpu_end-cpu_begin, cpu_end, cpu_begin
WRITE(Out_Unit,*)cpu_end,cpu_begin,ntheta, "pos5"

!---------------------------------------------------------------------
! Calculate chi DAF.
!---------------------------------------------------------------------
IF (debug) WRITE(Out_Unit, * ) 'T_chi'
IF (nchi>1)THEN
   CALL cputime (cpu_begin)
   Identifier="_chi"
   CALL get_daf (nchi, chi_val, dchi, mval, dafx, hermit, ndaf,  nd_daf, toe, usys2, kmax, Identifier)
   CALL cputime (cpu_end)
   WRITE(Out_Unit, * ) 'get_daf: time=', cpu_end-cpu_begin
ENDIF

!---------------------------------------------------------------------
! Calculate chi Kinetic energy term.
!---------------------------------------------------------------------
CALL cputime (cpu_begin)
CALL get_kinetic (nchi,   toe, t_chi  , chi_val, kmax, 'chi     ', system, group, tblock, iblock, usys2, dafx, nd_daf, ndaf)
CALL cputime (cpu_end)

WRITE(Out_Unit, * ) 'get_kinetic: time=', cpu_end-cpu_begin
!---------------------------------------------------------------------
! Calculate chi Coriolis term.
!---------------------------------------------------------------------
CALL cputime (cpu_begin)
!cc      CALL get_kinetic(nchi, toe, t_chi1, chi_val, kmax,'coriolis',
!cc     >                 system, group, tblock, iblock, usys2,dafx,
!cc     >                 nd_daf, ndaf)
CALL cputime (cpu_end)

WRITE(Out_Unit, * ) 'get_kinetic: time=', cpu_end-cpu_begin
CALL cputime (cpu_begin)
CALL symop (1, 1, 1, group, 0, nchi, 1, n1, nirrep, nsymop, 0, system, fac, par)
CALL cputime (cpu_end)
WRITE(Out_Unit, * ) 'symop: time=', cpu_end-cpu_begin
!---------------------------------------------------------------------
! Loop over each irreducible representation.
!---------------------------------------------------------------------
DO irrep = 1, nirrep
   CALL repnames (irrep, group, repname)
   WRITE(Out_Unit,  * ) 'irrep=', irrep, '  group=', group, '  representation=', repname, '  system=', system
   nchir = n1 (irrep)
   nsize = nrho * (ntheta) * nchir
   IF (irrep==1) npot = nsize
   nchan = min (nchanl, nsize)

   parity = par (irrep)
   parity = 0
   IF (jtot/=0.or.parity/=1)THEN
      !----------------------------------------------------------------------
      ! IF direct=.true. construct Hamiltonian and diagonalize.
      !----------------------------------------------------------------------
      IF (direct)THEN
         WRITE(Out_Unit, * ) 'Constructing Hamiltonian: nsize=', nsize
         WRITE(Out_Unit,  * ) 'irrep=', irrep, '   nirrep=', nirrep
         CALL cputime (cpu_begin)
         CALL hmlold (nrho, rho_val, ntheta, tht_val, nchi, nchir, &
          chi_val, t_rho, t_theta, t_chi, t_chir, t_chi1, t_chir1, &
          hamil, v_pot, nirrep, irrep, n1, nsize, jtot, parity, usys2, &
          a, b, c, d, rhom2, tblock)
         CALL cputime (cpu_end)
         WRITE(Out_Unit, * ) 'hml: time=', cpu_end-cpu_begin
         IF (debug)THEN
            WRITE(Out_Unit, * ) 'Hamiltonian Matrix'
            CALL MxOut(hamil, nsize, nsize)
         ENDIF
         fulldiag = .True.
         IF (fulldiag)THEN
            WRITE(Msg_Unit, * ) 'calling devcsf'
            CALL devcsf (nsize, hamil, nsize, eig, psi, nsize) !tempmod
            WRITE(Msg_Unit, * ) 'calling sorter'
            CALL sorter (eig, eig_im, nsize, psi)
            WRITE(Msg_Unit, * ) 'returned from sorter'
         ELSE
            WRITE(Msg_Unit, * ) 'calling devesf', nsize, nchanl
            !CALL devesf (nsize, min (nchanl, nsize), hamil, nsize, .true., eig, psi, nsize) !tempmod
            !STOP 'must supply routine'
            RETURN !TEMPMOD  
         ENDIF
         IF (debug) WRITE(Out_Unit, * ) 'info=', info
         !-------------------------------------------------------------------
         ! IF direct=.false. use iterative procedure to diagonalize
         ! Hamiltonian and never explicitly construct the Hamiltonian.
         !-------------------------------------------------------------------
      ELSE
         STOP 'indirect'
      ENDIF
      nbound = nchan
      WRITE(Out_Unit, * ) 'v_pot(1),v_pot(npot)', v_pot (1) , v_pot (npot)
      DO jval = 1, nchan
         IF (eig (jval) >v_pot (1) .or.eig (jval) >v_pot (npot) )THEN
            nbound = jval - 1
            goto 13
         ENDIF
      ENDDO
      !---------------------------------------------------------------------
      ! Print bound eigenvalues.
      !---------------------------------------------------------------------
      13    IF (irrep==1) eiga1 = eig (1)
      IF (nrho==1)THEN
         prefac = usys2 * rho_val (1) **2
         IF (ntheta==1)THEN
            prefac = usys2 * rho_val (1) **2 / (2 * b (1) )
         ENDIF
      ELSE
         prefac = usys2
      ENDIF
      WRITE(Out_Unit, * ) 'prefac=', prefac
      WRITE(Out_Unit, * ) 'There are ', nbound, ' bound eigenvalues(BoundOld)'
      IF (nbound/=0)THEN
         !WRITE(Out_Unit,*) 'Bound Eigenvalues in ev'
         WRITE(Out_Unit,*) 'Bound Eigenvalues in Hartrees'
         WRITE(Out_Unit,15) 'eig(jval)', 'eig(jval)-eiga1', 'prefac*eig(jval)', 'prefac*(eig(jval)-eiga1)'
         DO jval = 1, nbound
            !WRITE(Out_Unit,30) eig (jval)*autoev, (eig (jval) - eiga1)*autoev, prefac * eig (jval), prefac * (eig (jval) - eiga1)
            WRITE(Out_Unit,30) eig (jval), (eig (jval) - eiga1), prefac * eig (jval), prefac * (eig (jval) - eiga1)
         ENDDO
         WRITE(Out_Unit,*)
      ENDIF
      !---------------------------------------------------------------------
      ! Print scattering eigenvalues.
      !---------------------------------------------------------------------
      IF (nchan - nbound>0)THEN
         WRITE(Out_Unit,*) 'There are ', nchan - nbound, ' scattering eigenvalues'
         !WRITE(Out_Unit,*) 'Scattering Eigenvalues in ev'
         WRITE(Out_Unit,*) 'Scattering Eigenvalues in Hartrees'
         WRITE(Out_Unit,15) 'eig(jval)', 'eig(jval)-eiga1', 'prefac*eig(jval)', 'prefac*(eig(jval)-eiga1)'
         DO jval = nbound+1, nchan
            !WRITE(Out_Unit, 30) eig (jval)*autoev, (eig (jval) - eiga1)*autoev, prefac * eig (jval), prefac * (eig (jval) - eiga1)
            WRITE(Out_Unit, 30) eig (jval), (eig (jval) - eiga1), prefac * eig (jval), prefac * (eig (jval) - eiga1)
         ENDDO
         WRITE(Out_Unit,*)
      ENDIF
      15 FORMAT(24x,A9,11x,A15,7x,A16,4x,A24)
      !---------------------------------------------------------------------
      ! Store results for subsequent plotting.
      !---------------------------------------------------------------------
      lammax = jtot
      IF (parity==0)THEN
         lammin = 0
      ELSE
         lammin = 1
      ENDIF
      IF (irrep==1)THEN
         nbegin = 1
      ELSE
         nbegin = nbegin + n1 (irrep - 1)
      ENDIF
      nsize_new = nrho * (ntheta) * nchi * (lammax - lammin + 1)
      ALLOCATE (psi_new (nsize_new, nchan) )
      WRITE(Msg_Unit, * ) 'calling btrns'
      CALL btrns (nrho, ntheta, nchi, nchan, lammin, lammax, nchir, nsize, psi, nsize_new, psi_new, tblock (1, nbegin) )
      WRITE(Msg_Unit, * ) 'calling grafx in boundold'
      WRITE(Msg_Unit, * ) nrho, ntheta, nchi, nchan, nsize_new, repname
      CALL grafx(nrho, rho_val, ntheta, tht_val, nchi, chi_val, eig,   &
           v_pot, psi_new, nchan, nsize_new, repname)
      WRITE(Msg_Unit, * ) 'returned from calling grafx in boundold'
      DEALLOCATE (psi_new)
   ENDIF
ENDDO

WRITE(Out_Unit, * ) 'END(BoundOld)'

WRITE(Out_Unit, * )
RETURN
   30 FORMAT(1x,'Eigenvalues: ',4es23.12)  
ENDSUBROUTINE boundold
