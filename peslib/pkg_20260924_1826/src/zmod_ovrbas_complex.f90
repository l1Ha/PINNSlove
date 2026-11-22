module mod_ovrbas_complex
  implicit none
  private
  public :: ovrbas_complex

contains

  subroutine ovrbas_complex(rholast, rho, cvlast, cvnow, cprev, ccurr, ilev, jlev, klev, s, nprev, ncurr)
    use mod_common
    use mod_overlp, only: overlp
    implicit double precision (a-h,o-z)

    double precision, intent(in) :: rholast, rho
    integer, intent(in) :: nprev, ncurr
    double precision, intent(in) :: cvlast(nvi,n), cvnow(nvi,n)
    integer, intent(in) :: ilev(n), jlev(n), klev(n)
    complex*16, intent(in) :: cprev(nchan,nchan), ccurr(nchan,nchan)
    complex*16, intent(out) :: s(nchan,nchan)

    double precision, allocatable :: hover(:,:)
    complex*16, allocatable :: hcmp(:,:), tmp(:,:)
    integer :: i, j

    allocate(hover(nchan,nchan), hcmp(nchan,nchan), tmp(nchan,nchan))
    call overlp(rholast, cvlast, rho, cvnow, ilev, jlev, klev, hover)
    do j = 1, nchan
      do i = 1, nchan
        hcmp(i,j) = dcmplx(hover(i,j), 0.d0)
      enddo
    enddo
    call zgemm('N','N',nchan,ncurr,nchan,(1.d0,0.d0),hcmp,nchan,ccurr,nchan,(0.d0,0.d0),tmp,nchan)
    call zgemm('C','N',nprev,ncurr,nchan,(1.d0,0.d0),cprev,nchan,tmp,nchan,(0.d0,0.d0),s,nchan)
    deallocate(hover, hcmp, tmp)
  end subroutine ovrbas_complex

end module mod_ovrbas_complex
