module mod_sector_complex
  private
  public :: sector_complex

contains

  subroutine sector_complex(rho, cvi, cro, ilev, jlev, klev, nlev, h, s, c, e, np)
    use mod_sbasis, only: sbasis
    use mod_couple, only: couple
    use mod_hcsevp_complex, only: hcsevp_complex
    use mod_common
    implicit double precision (a-h,o-z)

    dimension cvi(nvi,n), cro(3,n)
    dimension ilev(n), jlev(n), klev(n), nlev(n)
    complex*16 h(nchan,nchan), s(nchan,nchan), c(nchan,nchan)
    double precision e(nchan)

    allocatable eint(:), sigma(:)
    allocatable hr(:,:), sr(:,:)
    integer :: i, j, ierr
    double precision, parameter :: cutoff = 1.d-4

    allocate(eint(n), sigma(nchan), hr(nchan,nchan), sr(nchan,nchan))

    call sbasis(rho, cvi, eint, ilev, jlev, klev, nlev, 0)
    call couple(rho, cvi, cro, eint, ilev, jlev, klev, hr, sr)

    do j = 1, nchan
      do i = 1, nchan
        h(i,j) = dcmplx(hr(i,j), 0.d0)
        s(i,j) = dcmplx(sr(i,j), 0.d0)
      enddo
    enddo

    call hcsevp_complex(h, s, c, nchan, e, sigma, np, cutoff, ierr)
    if (ierr .ne. 0) stop 'sector_complex 1'

    write(6,61) e(1)/rmlmda, e(np)/rmlmda, np
61  format(1x,'SECTORC:', '  e(1) = ',f8.3, '  e(np) = ',f8.3, '    np = ',i8)

    deallocate(eint, sigma, hr, sr)
  end subroutine sector_complex

end module mod_sector_complex
