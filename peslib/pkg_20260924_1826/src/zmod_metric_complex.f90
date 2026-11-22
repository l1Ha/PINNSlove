module mod_metric_complex
  private
  public :: metric_complex

contains

  subroutine metric_complex(rhop, cvp, npp, cprev, rho, cvi, np, ccurr, ilev, jlev, klev, h, s, t)
    use mod_common
    use mod_overlp, only: overlp
    use mod_scratchfile_complex, only: getrec_complex, putrec_complex
    implicit double precision (a-h,o-z)

    dimension cvp(nvi,n), cvi(nvi,n)
    dimension ilev(n), jlev(n), klev(n)
    complex*16 cprev(nchan,nchan), ccurr(nchan,nchan)
    complex*16 h(nchan,nchan), s(nchan,nchan), t(nchan,nchan)

    allocatable hover(:,:), htmp(:,:)
    integer :: i, j

    allocate(hover(nchan,nchan), htmp(nchan,nchan))
    call overlp(rhop, cvp, rho, cvi, ilev, jlev, klev, hover)

    do j = 1, nchan
      do i = 1, nchan
        h(i,j) = dcmplx(hover(i,j), 0.d0)
      enddo
    enddo

    call zgemm('N','N',nchan,np,nchan,(1.d0,0.d0),h,nchan,ccurr,nchan,(0.d0,0.d0),s,nchan)
    call getrec_complex(h, nchan, 0)
    call putrec_complex(ccurr, nchan, 0)
    call zgemm('C','N',npp,np,nchan,(1.d0,0.d0),cprev,nchan,s,nchan,(0.d0,0.d0),t,nchan)

    deallocate(hover, htmp)
  end subroutine metric_complex

end module mod_metric_complex
