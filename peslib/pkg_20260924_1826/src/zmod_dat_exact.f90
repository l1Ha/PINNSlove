module mod_dat_exact
  use mod_pes_multi, only: set_pes_lambda_scale
  use mod_ovrbas_complex, only: ovrbas_complex
  use mod_dat_overlap_chain, only: combine_dat_overlaps
  use mod_sector_complex, only: sector_complex
  implicit none
  private
  public :: dat_chain_exact

contains

  subroutine dat_chain_exact(rho, nsteps, cvi, cro, ilev, jlev, klev, nlev, overlaps, u, ierr)
    use mod_common
    implicit double precision (a-h,o-z)

    double precision, intent(in) :: rho
    integer, intent(in) :: nsteps
    double precision, intent(inout) :: cvi(nvi,n)
    double precision, intent(in) :: cro(3,n)
    integer, intent(inout) :: ilev(n), jlev(n), klev(n), nlev(n)
    complex*16, intent(out) :: overlaps(nchan,nchan,nsteps), u(nchan,nchan)
    integer, intent(out) :: ierr

    double precision :: lam_prev, lam_curr
    double precision, allocatable :: cv_prev(:,:), cv_curr(:,:)
    complex*16, allocatable :: c_prev(:,:), c_curr(:,:), h(:,:), s(:,:)
    double precision, allocatable :: e(:)
    integer :: np_prev, np_curr, istep

    allocate(cv_prev(nvi,n), cv_curr(nvi,n), c_prev(nchan,nchan), c_curr(nchan,nchan))
    allocate(h(nchan,nchan), s(nchan,nchan), e(nchan))

    ierr = 0
    overlaps = (0.d0, 0.d0)
    do istep = 1, nsteps
      do np_prev = 1, nchan
        overlaps(np_prev,np_prev,istep) = (1.d0, 0.d0)
      enddo
    enddo

    lam_prev = 1.d0
    call set_pes_lambda_scale(lam_prev)
    call sector_complex(rho, cvi, cro, ilev, jlev, klev, nlev, h, s, c_prev, e, np_prev)
    cv_prev = cvi

    do istep = 1, nsteps
      lam_curr = 1.d0 - dble(istep)/dble(nsteps)
      call set_pes_lambda_scale(lam_curr)
      call sector_complex(rho, cvi, cro, ilev, jlev, klev, nlev, h, s, c_curr, e, np_curr)
      cv_curr = cvi
      call ovrbas_complex(rho, rho, cv_prev, cv_curr, c_prev, c_curr, ilev, jlev, klev, overlaps(:,:,istep), np_prev, np_curr)
      c_prev = c_curr
      cv_prev = cv_curr
      np_prev = np_curr
    enddo

    ! The propagation after rdat should remain on the adiabatic surface (lambda = 0)
    ! so we do NOT reset pes_lambda_scale to 1.d0 here.
    
    call combine_dat_overlaps(overlaps, nchan, nsteps, u)
  end subroutine dat_chain_exact

end module mod_dat_exact
