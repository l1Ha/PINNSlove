module mod_solve_complex
  private
  public :: solve_complex

contains

  subroutine solve_complex(cvi, cro, ilev, jlev, klev, nlev, ered, nnrg)
    use mod_common, only: mro, mvi, nvi, n, nchan, rmin, rmax, mtr, rdat
    use mod_scratchfile_complex, only: putrec_complex, getrec_complex
    use mod_logder_complex, only: logder_complex
    use mod_sector_complex, only: sector_complex
    use mod_metric_complex, only: metric_complex
    use mod_dat_exact, only: dat_chain_exact
    use mod_dat_overlap_chain, only: transform_logder_complex
    use mod_sbasis, only: sbasis
    use mod_overlp, only: overlp
    use mod_pes_multi, only: NPES
    implicit double precision (a-h,o-z)

    dimension cvi(nvi,n), cro(3,n)
    dimension ilev(n), jlev(n), klev(n), nlev(n)
    dimension ered(nnrg)

    allocatable cvp(:,:)
    allocatable eint(:), e(:)
    complex*16, allocatable :: h(:,:), s(:,:), c(:,:), cprev(:,:)
    allocatable hreal(:,:)

    integer :: np, npp, ktr, j, k, nrg, ierr, nsteps_dat
    double precision :: drho, rhop, rho
    logical :: dat_done
    complex*16, allocatable :: overlaps(:,:,:), udat(:,:), ytmp(:,:)

    dat_done = .false.
    nsteps_dat = 8

    write(6,'(A,1x,I0,1x,A,1x,I0,1x,A,1x,I0)') 'SOLVE_COMPLEX: start mtr', mtr, 'nchan', nchan, 'nnrg', nnrg
    flush(6)

    allocate(cvp(nvi,n), eint(n), e(nchan))
    allocate(h(nchan,nchan), s(nchan,nchan), c(nchan,nchan), cprev(nchan,nchan))
    allocate(hreal(nchan,nchan))
    allocate(ytmp(nchan,nchan), overlaps(nchan,nchan,nsteps_dat), udat(nchan,nchan))

    npp = 0
    drho = (rmax-rmin)/mtr
    rhop = rmin - 0.5d0*drho
    do ktr = 1, mtr
      if (ktr .eq. 1 .or. mod(ktr,10) .eq. 0 .or. ktr .eq. mtr) then
        write(6,'(A,1x,I0,1x,A,1x,I0,1x,A,1x,ES12.4)') 'SOLVE_COMPLEX: sector', ktr, 'of', mtr, 'rho', rhop + drho
        flush(6)
      endif
      rho = rhop + drho
      call sector_complex(rho, cvi, cro, ilev, jlev, klev, nlev, h, s, c, e, np)
      if (npp .eq. 0) then
        call putrec_complex(c, nchan, 0)
      else
        call metric_complex(rhop, cvp, npp, cprev, rho, cvi, np, c, ilev, jlev, klev, h, s, c)
      endif
      call logder_complex(drho, ered, nnrg, h, s, c, e, nchan, np, npp)

      if (.not. dat_done .and. NPES .gt. 1 .and. rdat .gt. 0.d0) then
        if (.not. dat_done .and. rho .ge. rdat) then
          ! --- Diabatic-to-Adiabatic (DAT) Transformation ---
          ! This step uncouples the electronic states at hyperradius rdat.
          ! Based on the theory described in [Kendrick, JCP 142, 054110 (2015)].
          ! See Section 2 (Methodology) of the accompanying article.
          write(6,'(A,1x,ES12.4,1x,A,1x,I0)') 'SOLVE_COMPLEX: exact DAT transform at rdat', rdat, 'steps', nsteps_dat
          flush(6)
          call dat_chain_exact(rdat, nsteps_dat, cvi, cro, ilev, jlev, klev, nlev, overlaps, udat, ierr)
          if (ierr .ne. 0) stop 'solve_complex 1'

          ! Transform the log-derivative matrix to the adiabatic basis
          do nrg = 1, nnrg
            call getrec_complex(ytmp, nchan, nrg)
            call transform_logder_complex(ytmp, udat, nchan, ytmp)
            call putrec_complex(ytmp, nchan, nrg)
          enddo
          call transform_logder_complex(cprev, udat, nchan, cprev)
          dat_done = .true.
          write(6,'(A)') 'SOLVE_COMPLEX: exact DAT transform done'
          flush(6)
        endif
      endif
      do j = 1, n
        do k = 1, nvi
          cvp(k,j) = cvi(k,j)
        enddo
      enddo
      cprev = c
      npp = np
      rhop = rho
    enddo

    npp = nchan + 1
    rho = rmax
    write(6,'(A,1x,ES12.4)') 'SOLVE_COMPLEX: final asymptotic transform at rho', rho
    flush(6)
    call sbasis(rho, cvi, eint, ilev, jlev, klev, nlev, 1)
    call overlp(rhop, cvp, rho, cvi, ilev, jlev, klev, hreal)
    do j = 1, nchan
      do k = 1, nchan
        h(k,j) = dcmplx(hreal(k,j), 0.d0)
      enddo
    enddo
    call getrec_complex(s, nchan, 0)
    call zgemm('C','N',nchan,np,nchan,(1.d0,0.d0),h,nchan,s,nchan,(0.d0,0.d0),c,nchan)
    call logder_complex(drho, ered, nnrg, h, s, c, e, nchan, np, npp)
    write(6,'(A)') 'SOLVE_COMPLEX: done'
    flush(6)
  end subroutine solve_complex

end module mod_solve_complex
