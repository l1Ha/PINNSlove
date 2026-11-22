module mod_pes_multi
  implicit double precision (a-h,o-z)
  public :: NPES, PESMODE, init_pes, potenl_matrix, potenl_ref
  public :: active_npes, asymptotic_npes
  public :: set_pes_lambda_scale, get_pes_lambda_scale

  integer :: NPES = 1
  character(len=20) :: PESMODE = 'INTERNAL_1STATE'
  real(8), parameter :: AU_TO_EV = 27.211386245988d0
  real(8) :: pes_lambda_scale = 1.d0

contains

  subroutine init_pes(npes_in, mode)
    use mod_common, only: pes_C6
    integer, intent(in) :: npes_in
    character(len=*), intent(in) :: mode
    NPES = npes_in
    PESMODE = adjustl(mode)
    if (trim(PESMODE) .eq. 'KRb2' .or. trim(PESMODE) .eq. 'KRB2') then
      call krb2_pes_init
      if (maxval(abs(pes_C6)) .le. 0.d0) then
        pes_C6 = (/ 4299.509938d0, 4299.509938d0, 4710.216338d0 /)
      end if
    else if (trim(PESMODE) .eq. 'H3') then
      write(6,*) 'init_pes: H3 PES mode, NPES=', NPES
    end if
  end subroutine init_pes

  subroutine set_pes_lambda_scale(scale)
    real(8), intent(in) :: scale
    pes_lambda_scale = scale
  end subroutine set_pes_lambda_scale

  real(8) function get_pes_lambda_scale()
    get_pes_lambda_scale = pes_lambda_scale
  end function get_pes_lambda_scale

  integer function active_npes(rho)
    real(8), intent(in) :: rho

    active_npes = max(1, NPES)
  end function active_npes

  integer function asymptotic_npes()
    asymptotic_npes = max(1, NPES)
  end function asymptotic_npes

  real(8) function pair_distance(i, j, rvec)
    integer, intent(in) :: i, j
    real(8), intent(in) :: rvec(3)

    if ((i .eq. 1 .and. j .eq. 2) .or. (i .eq. 2 .and. j .eq. 1)) then
      pair_distance = rvec(3)
    else if ((i .eq. 1 .and. j .eq. 3) .or. (i .eq. 3 .and. j .eq. 1)) then
      pair_distance = rvec(2)
    else if ((i .eq. 2 .and. j .eq. 3) .or. (i .eq. 3 .and. j .eq. 2)) then
      pair_distance = rvec(1)
    else
      pair_distance = 0.d0
    end if
  end function pair_distance

  subroutine potenl_matrix(ra, sa, cosa, V, ia)
    use mod_common
    use mod_mkmp2, only: pot0
    implicit double precision (a-h,o-z)
    real(8), intent(in) :: ra, sa, cosa
    integer, intent(in) :: ia
    real(8), intent(out) :: V(:,:)
    integer :: i, ib, ic, ik, ir1, ir2
    real(8) :: R, E1inf, E2inf, D1, D2, R01, R02, sigma1, sigma2
    real(8) :: Vc, Rc, sigmac, tmp1, tmp2
    real(8) :: rvec(3), q(3), vd11, vd22, vd12, Va1, Va2, repulsive_wall
    real(8) :: r1, r2, r3
    real(8) :: rap, cap, sap, sbp, scp
    real(8) :: vdia(2,2)

    V = 0.d0

    if (trim(PESMODE) .eq. 'H3') then
      ib = ia + 1
      if (ib .gt. 3) ib = 1
      ic = 6 - ia - ib

      rap = ra/scales(ia)
      cap = 2.0d0*cosa*rap
      sap = scales(ia)*sa
      sbp = mass(ib)/(mass(ib)+mass(ic))*sap
      scp = mass(ic)/(mass(ib)+mass(ic))*sap

      rvec(ia) = sap
      rvec(ib) = sqrt(max(1.d-10, rap**2 - cap*sbp + sbp**2))
      rvec(ic) = sqrt(max(1.d-10, rap**2 + cap*scp + scp**2))

      if (rvec(1).le.0.5d0 .or. rvec(2).le.0.5d0 .or. rvec(3).le.0.5d0) then
        V(1,1) = AU_TO_EV * 3.d0
        if (NPES .ge. 2) then
          V(2,2) = AU_TO_EV * 3.d0
          V(1,2) = 0.d0
          V(2,1) = 0.d0
        endif
        return
      endif

      call pot0(NPES, vdia, rvec, pes_lambda_scale)
      V(1,1) = AU_TO_EV * vdia(1,1)
      if (NPES .ge. 2) then
        V(2,2) = AU_TO_EV * vdia(2,2)
        V(1,2) = AU_TO_EV * vdia(1,2)
        V(2,1) = AU_TO_EV * vdia(2,1)
      end if
      return

    else if (NPES .eq. 1 .and. (trim(PESMODE) .eq. 'KRb2' .or. trim(PESMODE) .eq. 'KRB2')) then
      ib = ia + 1
      if (ib .gt. 3) ib = 1
      ic = 6 - ia - ib

      rap = ra/scales(ia)
      cap = 2.0d0*cosa*rap
      sap = scales(ia)*sa
      sbp = mass(ib)/(mass(ib)+mass(ic))*sap
      scp = mass(ic)/(mass(ib)+mass(ic))*sap

      rvec(ia) = sap
      rvec(ib) = sqrt(max(0.d0, rap**2 - cap*sbp + sbp**2))
      rvec(ic) = sqrt(max(0.d0, rap**2 + cap*scp + scp**2))

      ik = 1
      do i = 2,3
        if (mass(i) .lt. mass(ik)) ik = i
      enddo
      ir1 = 1
      if (ir1 .eq. ik) ir1 = 2
      ir2 = 6 - ik - ir1
      q(1) = pair_distance(ik, ir1, rvec)
      q(2) = pair_distance(ik, ir2, rvec)
      q(3) = pair_distance(ir1, ir2, rvec)

      call krb2_pes_eval(q, pes_lambda_scale, vd11, vd22, vd12, Va1, Va2)
      V(1,1) = AU_TO_EV * Va1
      return
    else if (NPES .eq. 1) then
      V(1,1) = 0.0d0
      return
    else if (NPES .eq. 2 .and. trim(PESMODE) .eq. 'INTERNAL_2STATE') then
      R = ra

      ! default 2-state diabatic model parameters (can be changed later)
      E1inf = 0.0d0
      E2inf = 0.20d0
      D1 = -0.30d0
      D2 = -0.10d0
      R01 = 4.0d0
      R02 = 6.0d0
      sigma1 = 1.0d0
      sigma2 = 1.0d0
      Vc = 0.05d0
      Rc = 5.0d0
      sigmac = 0.8d0

      tmp1 = exp(-0.5d0*((R-R01)/sigma1)**2)
      tmp2 = exp(-0.5d0*((R-R02)/sigma2)**2)

      V(1,1) = E1inf + D1*tmp1
      V(2,2) = E2inf + D2*tmp2
      V(1,2) = Vc*exp(-0.5d0*((R-Rc)/sigmac)**2)
      V(2,1) = V(1,2)
      return
    else if (NPES .eq. 2 .and. (trim(PESMODE) .eq. 'KRb2' .or. trim(PESMODE) .eq. 'KRB2')) then
      ! Mapping of mass-scaled Jacobi coordinates (ra, sa, cosa) to interatomic distances.
      ! Based on the kinematic scaling theory by F. T. Smith, J. Chem. Phys. 31, 1352 (1959).
      ib = ia + 1
      if (ib .gt. 3) ib = 1
      ic = 6 - ia - ib

      rap = ra/scales(ia)  ! Real-space Jacobi R
      cap = 2.0d0*cosa*rap
      sap = scales(ia)*sa  ! Real-space Jacobi r
      sbp = mass(ib)/(mass(ib)+mass(ic))*sap
      scp = mass(ic)/(mass(ib)+mass(ic))*sap

      rvec(ia) = sap
      rvec(ib) = sqrt(max(0.1d0, rap**2 - cap*sbp + sbp**2))
      rvec(ic) = sqrt(max(0.1d0, rap**2 + cap*scp + scp**2))
      
      ! Interatomic distances q(1:3) following KRb2 PES convention (KRb, KRb, RbRb)
      ! Numerical protection: Cap minimum distance to 0.5 bohr to avoid PES singularities.
      ik = 1
      do i = 2,3
        if (mass(i) .lt. mass(ik)) ik = i
      enddo
      ir1 = 1
      if (ir1 .eq. ik) ir1 = 2
      ir2 = 6 - ik - ir1
      
      ! Prevent singularity at short range (ra -> 0)
      q(1) = max(0.5d0, pair_distance(ik, ir1, rvec))
      q(2) = max(0.5d0, pair_distance(ik, ir2, rvec))
      q(3) = max(0.5d0, pair_distance(ir1, ir2, rvec))

      call krb2_pes_eval(q, 1.0d0, vd11, vd22, vd12, Va1, Va2)
      vd11 = AU_TO_EV * vd11
      vd22 = AU_TO_EV * vd22
      vd12 = AU_TO_EV * vd12
      Va1 = AU_TO_EV * Va1
      Va2 = AU_TO_EV * Va2

      ! Add an artificial repulsive wall to prevent RKHS negative divergence at short range
      repulsive_wall = 0.0d0
      do i = 1, 3
         if (q(i) .lt. 5.5d0) then
            ! A smooth C12 repulsive core that ramps up strongly for q < 5.5 bohr
            repulsive_wall = repulsive_wall + 10.0d0 * ((5.5d0 - q(i)) / 1.0d0)**4
         endif
      enddo
      vd11 = vd11 + repulsive_wall
      vd22 = vd22 + repulsive_wall
      Va1 = Va1 + repulsive_wall
      Va2 = Va2 + repulsive_wall

      ! Smoothly interpolate between Adiabatic (lambda=0) and Diabatic (lambda=1)
      V(1,1) = pes_lambda_scale * vd11 + (1.0d0 - pes_lambda_scale) * Va1
      V(2,2) = pes_lambda_scale * vd22 + (1.0d0 - pes_lambda_scale) * Va2
      V(1,2) = pes_lambda_scale * vd12
      V(2,1) = V(1,2)
      
      ! Hard cap to avoid any NaN or extreme negative values
      if (V(1,1) .lt. -50.0d0) V(1,1) = 500.0d0
      if (V(2,2) .lt. -50.0d0) V(2,2) = 500.0d0
      return
    else
      do i = 1, size(V,1)
        V(i, :) = 0.0d0
      end do
      return
    end if
  end subroutine potenl_matrix

  subroutine potenl_ref(ra, sa, cosa, vref, ia)
    implicit double precision (a-h,o-z)
    real(8), intent(in) :: ra, sa, cosa
    integer, intent(in) :: ia
    real(8), intent(out) :: vref
    integer :: n
    real(8), allocatable :: V(:,:)
    real(8) :: a,b,c,disc,lam1,lam2

    n = max(1, NPES)
    allocate(V(n,n))
    call potenl_matrix(ra, sa, cosa, V, ia)

    if (NPES .eq. 1) then
      vref = V(1,1)
    else if (NPES .eq. 2) then
      a = V(1,1)
      b = V(1,2)
      c = V(2,2)
      disc = sqrt(max(0.0d0, (a-c)**2 + 4.0d0*b*b))
      lam1 = 0.5d0*(a + c - disc)
      lam2 = 0.5d0*(a + c + disc)
      vref = min(lam1, lam2)
    else
      vref = V(1,1)
    end if

    deallocate(V)
  end subroutine potenl_ref

  real(8) function get_elec_offset(ia, pp)
    implicit none
    integer, intent(in) :: ia, pp
    real(8) :: V(2,2), ra_inf
    
    get_elec_offset = 0.d0
    if (NPES .eq. 2 .and. pp .eq. 2) then
       ra_inf = 50.0d0
       call potenl_matrix(ra_inf, 1.0d0, 0.d0, V, ia)
       get_elec_offset = V(2,2) - V(1,1)
       if (get_elec_offset .lt. 0.d0) get_elec_offset = 1.0d5
    endif
  end function get_elec_offset

end module mod_pes_multi
