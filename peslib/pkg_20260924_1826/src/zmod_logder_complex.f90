module mod_logder_complex
  use mod_complex_linalg, only: solve_linear_complex
  use mod_scratchfile_complex, only: getrec_complex, putrec_complex
  implicit none
  private
  public :: logder_complex

contains

  subroutine invert_complex(a, n, ierr)
    integer, intent(in) :: n
    complex*16, intent(inout) :: a(n,n)
    integer, intent(out) :: ierr
    complex*16, allocatable :: rhs(:,:), amat(:,:)
    integer :: i

    allocate(rhs(n,n), amat(n,n))
    rhs = (0.d0,0.d0)
    do i = 1, n
       rhs(i,i) = (1.d0,0.d0)
    enddo
    amat = a
    call solve_linear_complex(amat, rhs, n, n, ierr)
    a = rhs
    deallocate(rhs, amat)
  end subroutine invert_complex

  subroutine logder_complex(drho, ered, nnrg, y, z, t, e, n, np, npp)
    ! Implements the multichannel log-derivative method for quantum scattering.
    ! Based on B. R. Johnson, J. Comput. Phys. 13, 445 (1973).
    integer, intent(in) :: nnrg, n, np, npp
    double precision, intent(in) :: drho, ered(nnrg), e(n)
    complex*16, intent(inout) :: y(n,n), z(n,n), t(n,n)

    integer :: nrg, i, j, ierr, npuse, nppuse
    double precision :: h, eps, psq, p
    double precision :: y1(n), y2(n), y3(n), y4(n)
    complex*16, allocatable :: yloc(:,:), zloc(:,:), tloc(:,:), tmp(:,:), inva(:,:)

    npuse = np
    nppuse = npp
    do nrg = 1, nnrg
       if (nppuse .le. n) then
          h = drho
          eps = 1.d-16/(h*h)
          do i = 1, npuse
             psq = e(i) - ered(nrg)
             if (psq .gt. eps) then
                p = sqrt(psq)
                y1(i) = p/tanh(p*h)
                y2(i) = p/sinh(p*h)
             else if (psq .lt. -eps) then
                p = sqrt(-psq)
                y1(i) = p/tan(p*h)
                y2(i) = p/sin(p*h)
             else
                y1(i) = 1.d0/h
                y2(i) = 1.d0/h
             endif
             y3(i) = y2(i)
             y4(i) = y1(i)
          enddo
          if (npuse .lt. n) then
             do i = npuse+1, n
                y1(i) = 0.d0
                y2(i) = 0.d0
                y3(i) = 0.d0
                y4(i) = 0.d0
             enddo
          endif
       endif

       if (nppuse .le. 0) then
          y = (0.d0,0.d0)
          do j = 1, n
             y(j,j) = dcmplx(y4(j), 0.d0)
          enddo
       else if (nppuse .le. n) then
          call getrec_complex(y, n, nrg)
          allocate(yloc(npuse,npuse), zloc(nppuse,npuse), tloc(nppuse,npuse), tmp(nppuse,npuse))
          zloc = (0.d0,0.d0)
          tloc = t(1:nppuse,1:npuse)
          call zgemm('N','N',nppuse,npuse,nppuse,(1.d0,0.d0),y,n,t,n,(0.d0,0.d0),zloc,nppuse)
          yloc = (0.d0,0.d0)
          call zgemm('C','N',npuse,npuse,nppuse,(1.d0,0.d0),t,n,zloc,nppuse,(0.d0,0.d0),yloc,npuse)
          do j = 1, npuse
             yloc(j,j) = yloc(j,j) + dcmplx(y1(j),0.d0)
          enddo
          call invert_complex(yloc, npuse, ierr)
          if (ierr .ne. 0) stop 'logder_complex 1'
          y = (0.d0,0.d0)
          do j = 1, npuse
             do i = 1, npuse
                y(i,j) = -dcmplx(y3(i),0.d0) * yloc(i,j) * dcmplx(y2(j),0.d0)
             enddo
             y(j,j) = y(j,j) + dcmplx(y4(j),0.d0)
          enddo
          deallocate(yloc, zloc, tloc, tmp)
       else
          call getrec_complex(y, n, nrg)
          allocate(zloc(npuse,n), tloc(npuse,n), yloc(n,n))
          tloc = t(1:npuse,1:n)
          zloc = (0.d0,0.d0)
          call zgemm('N','C',npuse,n,npuse,(1.d0,0.d0),y,n,t,n,(0.d0,0.d0),zloc,npuse)
          yloc = (0.d0,0.d0)
          call zgemm('N','N',n,n,npuse,(1.d0,0.d0),t,n,zloc,npuse,(0.d0,0.d0),yloc,n)
          y = yloc
          deallocate(zloc, tloc, yloc)
       endif

       call putrec_complex(y, n, nrg)
    enddo
  end subroutine logder_complex

end module mod_logder_complex
