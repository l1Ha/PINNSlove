module mod_prot_complex
  private
  public :: prot_complex

contains

  subroutine prot_complex(p, jmax, kmax, x, phi)
    implicit double precision (a-h,o-z)
    integer, intent(in) :: jmax, kmax
    double precision, intent(in) :: x, phi
    complex*16, intent(out) :: p(0:jmax,0:kmax)

    double precision :: pr(0:jmax,0:kmax)
    integer :: j, k

    call prot_real(pr, jmax, kmax, x)
    do k = 0, kmax
      do j = 0, jmax
        p(j,k) = dcmplx(pr(j,k), 0.d0) * exp(dcmplx(0.d0, dble(k)*phi))
      enddo
    enddo
  contains

    subroutine prot_real(p, jmax, kmax, x)
      implicit double precision (a-h,o-z)
      integer, intent(in) :: jmax, kmax
      double precision, intent(in) :: x
      double precision, intent(out) :: p(0:jmax,0:kmax)
      integer, parameter :: jtop = 100
      double precision, save :: a(0:jtop,0:jtop), b(0:jtop,0:jtop), c(0:jtop,0:jtop)
      integer, save :: init = 0
      integer :: j, k

      if (init .eq. 0) then
        do j = 0, jtop
          do k = 0, j-2
            a(j,k) = (j+j-1.d0)/(j-k)
            b(j,k) = (j+k-1.d0)/(j-k)
          enddo
          fac = j+0.5d0
          c(j,0) = sqrt(fac)
          do k = 1, j
            fac = fac/((j+k)*(j-k+1))
            c(j,k) = sqrt(fac)
          enddo
        enddo
        init = 1
      endif

      if (jmax .gt. jtop) stop 'prot_complex 0'
      if (abs(x) .gt. 1.d0) stop 'prot_complex 1'

      do j = 0, jmax
        do k = j+1, kmax
          p(j,k) = 0.d0
        enddo
      enddo
      p(0,0) = 1.d0
      if (kmax .gt. 0) then
        sx = -sqrt((1.d0-x)*(1.d0+x))
        do k = 1, min(jmax,kmax)
          p(k,k) = (2*k-1)*sx*p(k-1,k-1)
        enddo
      endif
      if (jmax .gt. 0) then
        p(1,0) = x
        do k = 1, min(jmax-1,kmax)
          p(k+1,k) = (2*k+1)*x*p(k,k)
        enddo
      endif
      do k = 0, kmax
        do j = k+2, jmax
          p(j,k) = a(j,k)*x*p(j-1,k)-b(j,k)*p(j-2,k)
        enddo
      enddo
      do k = 0, kmax
        do j = k, jmax
          p(j,k) = c(j,k)*p(j,k)
        enddo
      enddo
    end subroutine prot_real

  end subroutine prot_complex

end module mod_prot_complex
