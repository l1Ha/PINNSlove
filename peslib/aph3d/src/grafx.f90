SUBROUTINE grafx (nrho, rho_val, ntheta, theta_val, nchi, chi_val, eig, v_pot, psi, nchanl, npsi, repname)
USE Numeric_Kinds_Module
USE FileUnits_Module
USE DiatomicPot_Module, ONLY: DiatomicPot
!-----------------------------------------------------------------------
! This routine was written by G. A. Parker
! IF you find an error or have an improvement please send a messge to
! Parker@ou.edu
!-----------------------------------------------------------------------
IMPLICIT NONE
CHARACTER(LEN=3) :: repname, squish * 10
CHARACTER(LEN=40) :: fn1, fn2, fn3, flname * 17
INTEGER :: p1, p2, p3, ichanl
INTEGER :: nrho, ntheta, nchi, nchanl, npsi, ntp (100)
INTEGER, ALLOCATABLE::poly (:), poly_p (:)
REAL(Kind=WP_Kind) :: rho_val (nrho), theta_val (ntheta), chi_val (nchi)
REAL(Kind=WP_Kind) :: eig (nchanl), v_pot (nrho, ntheta, nchi)
REAL(Kind=WP_Kind) :: psi (npsi, nchanl)
REAL(KIND=WP_Kind), ALLOCATABLE::vert (:, :), vert_p (:, :), fx (:), fy (:),f (:)

WRITE(Out_Unit, * ) 'nrho,ntheta,nchi,nchanl=', nrho, ntheta, nchi, nchanl
WRITE(Out_Unit, * ) 'npsi=', npsi

WRITE(Out_Unit, * ) 'repname=', repname
ALLOCATE (vert (3, nrho * ntheta * nchi), vert_p (3, nchi) )
ALLOCATE (fx (nrho * ntheta * nchi), fy (nrho * ntheta * nchi) )
ALLOCATE (f (nrho * ntheta * nchi) )
ALLOCATE (poly (10 * nrho * ntheta * nchi), poly_p (nrho * nchi + 1) )
!-----------------------------------------------------------------------
! store the calculated wavefunctions for plotting.
! calculate turning points and store for plotting.
!-----------------------------------------------------------------------
IF (ntheta==1.AND.nchi==1)THEN
   fn1 = 'DiatomicPlots/'//TRIM(DiatomicPot)//'_rho_psif_graph.csv'
   fn2 = 'DiatomicPlots/'//TRIM(DiatomicPot)//'_rho_epsi_graph.csv'
   fn3 = 'DiatomicPlots/'//TRIM(DiatomicPot)//'_rho_tpts_graph.csv'
   CALL psi_graph (rho_val, psi, nrho, nchanl, fn1)
   CALL epsi_graph (rho_val, psi, nrho, nchanl, eig, fn2)
   CALL tp_graph (nchanl, nrho, rho_val, v_pot, eig, ntp, fn3)
ENDIF

IF (nrho==1.AND.nchi==1)THEN
   fn1 = 'DiatomicPlots/'//TRIM(DiatomicPot)//'_tht_psif_graph.csv'
   fn2 = 'DiatomicPlots/'//TRIM(DiatomicPot)//'_tht_epsi_graph.csv'
   fn3 = 'DiatomicPlots/'//TRIM(DiatomicPot)//'_tht_tpts_graph.csv'
   CALL psi_graph (theta_val, psi, ntheta, nchanl, fn1)
   CALL epsi_graph (theta_val, psi, ntheta, nchanl, eig, fn2)
   CALL tp_graph (nchanl, ntheta, theta_val, v_pot, eig, ntp, fn3)
ENDIF

IF (nrho==1.AND.ntheta==1)THEN
   WRITE(Out_Unit, * ) repname
   fn1 = 'DiatomicPlots/'//TRIM(DiatomicPot)//repname//'_psif_graph.csv'
   fn2 = 'DiatomicPlots/'//TRIM(DiatomicPot)//repname//'_epsi_graph.csv'
   fn3 = 'DiatomicPlots/'//TRIM(DiatomicPot)//repname//'_tpts_graph.csv'
   CALL psi_graph (chi_val, psi, nchi, nchanl, fn1)
   CALL epsi_graph (chi_val, psi, nchi, nchanl, eig, fn2)
   CALL tp_graph (nchanl, nchi, chi_val, v_pot, eig, ntp, fn3)
ENDIF

!-----------------------------------------------------------------------
! IF two dimensional problem store wavefunctions and potential for plott
!-----------------------------------------------------------------------
IF (nrho>1.AND.nchi>1.AND.ntheta==1)THEN
   squish = '0123456789'
   flname (1:15) = 'vrho_chi000.obj'
   p3 = 8 - theta_val (1) / (3.1415d0 / ntheta)
   flname (11:11) = squish (p3:p3)
   WRITE(Msg_Unit,  * ) flname, '   location 1 -- potential'
   CALL rho_chi (nrho, ntheta, nchi, rho_val, chi_val, v_pot, flname, vert, poly)
ENDIF

IF (nrho>1.AND.ntheta>1.AND.nchi==1)THEN
   squish = '0123456789'
   flname (1:17) = 'vrho_theta000.obj'
   p3 = 1 + theta_val (1) / (3.1415d0 / 8)
   flname (13:13) = squish (p3:p3)
   CALL rho_theta (nrho, ntheta, nchi, rho_val, theta_val, v_pot, flname, vert, poly)
ENDIF

IF (ntheta>1.AND.nchi>1.AND.nrho==1)THEN
   p3 = 1
   p2 = 1
   p1 = rho_val (1) + 1
   IF (p1>10)THEN
      p1 = p1 - 10
   ENDIF

   DO ichanl = 1, nchanl
   flname (2:10) = 'theta_chi'
   squish = '0123456789'
   flname (11:11) = squish (p1:p1)
   flname (12:12) = squish (p2:p2)
   flname (13:13) = squish (p3:p3)
   flname (14:17) = '.obj'
   flname (1:1) = 'p'
   WRITE(Msg_Unit,  * ) flname, '   location 2 -- wavefunction'

   CALL theta_chi (nrho, ntheta, nchi, theta_val, chi_val, psi (1, &
    ichanl), flname, eig, ichanl, nchanl, vert, poly, vert_p,  poly_p, fx, fy, f)

   IF (ichanl==1)THEN
      flname (1:1) = 'v'
      WRITE(Msg_Unit,  * ) flname, '   location 3 -- potential'
      CALL theta_chi (nrho, ntheta, nchi, theta_val, chi_val, &
       v_pot, flname, eig, ichanl, nchanl, vert, poly, vert_p,  poly_p, fx, fy, f)
   ENDIF
   
   p3 = p3 + 1
   IF (p3==11)THEN
      p3 = 1
      p2 = p2 + 1
   ENDIF
   
   IF (p2==11)THEN
      p2 = 1
      p1 = p1 + 1

   ENDIF
   
   ENDDO
ENDIF

DEALLOCATE (vert, vert_p)
DEALLOCATE (fx, fy, f)

DEALLOCATE (poly, poly_p)
RETURN
ENDSUBROUTINE grafx
