module mod_hcsevp_complex
  use mod_complex_linalg, only: hermitian_eig_complex
  implicit none
  private
  public :: hcsevp_complex

contains

  subroutine hcsevp_complex(h, s, c, n, e, sigma, np, cutoff, ierr)
    integer, intent(in) :: n
    complex*16, intent(inout) :: h(n,n), s(n,n)
    complex*16, intent(out) :: c(n,n)
    double precision, intent(out) :: e(n), sigma(n)
    double precision, intent(in) :: cutoff
    integer, intent(out) :: np, ierr

    integer :: i, j, idx
    complex*16, allocatable :: t(:,:), horth(:,:)
    integer, allocatable :: keep(:)

    call hermitian_eig_complex(s, n, sigma, ierr)
    if (ierr .ne. 0) return

    allocate(keep(n))
    np = 0
    do j = 1, n
       if (sigma(j) .gt. cutoff) then
          np = np + 1
          keep(np) = j
       endif
    enddo
    if (np .le. 0) then
        deallocate(keep)
        ierr = -1
        return
    endif

    allocate(t(n,np), horth(np,np))
    do j = 1, np
       idx = keep(j)
       sigma(j) = sigma(idx)
       do i = 1, n
          t(i,j) = s(i,idx) / sqrt(sigma(j))
       enddo
    enddo

    call zgemm('C','N',np,n,np,(1.d0,0.d0),t,n,h,n,(0.d0,0.d0),c,n)
    call zgemm('N','N',np,np,n,(1.d0,0.d0),c,n,t,n,(0.d0,0.d0),horth,np)
    call hermitian_eig_complex(horth, np, e, ierr)
    if (ierr .ne. 0) then
       deallocate(t, horth, keep)
       return
    endif

    c = (0.d0, 0.d0)
    call zgemm('N','N',n,np,np,(1.d0,0.d0),t,n,horth,np,(0.d0,0.d0),c,n)
    deallocate(t, horth, keep)
  end subroutine hcsevp_complex

end module mod_hcsevp_complex
