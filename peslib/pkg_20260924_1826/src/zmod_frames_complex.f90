module mod_frames_complex
  use mod_frames, only: frames
  implicit none
  private
  public :: frames_complex

contains

  subroutine frames_complex(x, c, cro, jlev, klev, llev, mode)
    use mod_common
    implicit double precision (a-h,o-z)

    complex*16, intent(inout) :: x(nchan,nchan)
    double precision, intent(out) :: c(nchan,nchan), llev(nchan)
    double precision, intent(in) :: cro(3,n)
    integer, intent(in) :: jlev(n), klev(n), mode
    double precision, allocatable :: xr(:,:), xi(:,:)
    integer :: i, j

    allocate(xr(nchan,nchan), xi(nchan,nchan))
    do j = 1, nchan
      do i = 1, nchan
        xr(i,j) = dble(x(i,j))
        xi(i,j) = dimag(x(i,j))
      enddo
    enddo
    call frames(xr, c, cro, jlev, klev, llev, mode)
    call frames(xi, c, cro, jlev, klev, llev, mode)
    do j = 1, nchan
      do i = 1, nchan
        x(i,j) = dcmplx(xr(i,j), xi(i,j))
      enddo
    enddo
    deallocate(xr, xi)
  end subroutine frames_complex

end module mod_frames_complex
