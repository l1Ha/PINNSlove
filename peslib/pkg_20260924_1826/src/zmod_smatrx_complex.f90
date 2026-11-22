module mod_smatrx_complex
  use mod_asymptotic_driver_complex, only: asymptotic_driver_complex
  use mod_complex_linalg, only: hermitian_eig_complex
  implicit none
  private
  public :: smatrx_complex

contains

  subroutine smatrx_complex(kmat, dkde, ederiv, smat, qmat, qeval, eint, llev, n, seps, ierr, ered, phase_full, nopen_out)
    implicit double precision (a-h,o-z)
    integer, intent(in) :: n
    complex*16, intent(in) :: kmat(n,n)
    complex*16, intent(in) :: dkde(n,n)
    logical, intent(in) :: ederiv
    complex*16, intent(out) :: smat(n,n), qmat(n,n)
    double precision, intent(in) :: eint(n), llev(n), ered
    double precision, intent(out) :: seps, qeval(n), phase_full(n)
    integer, intent(out) :: ierr
    integer, intent(out) :: nopen_out

    complex*16, allocatable :: kop(:,:), sop(:,:), qop(:,:), dk(:,:), tran(:,:)
    complex*16, allocatable :: qwork(:,:)
    double precision, allocatable :: phase(:)
    integer :: i, j, ic, jc, nopen
    double precision :: sumv, difv, pi
    complex*16 :: cphase, qphase

    sumv = 1.d0
    difv = 0.d0
    nopen = 0
    pi = acos(-1.d0)
    smat = (0.d0,0.d0)
    qmat = (0.d0,0.d0)
    qeval = 0.d0
    phase_full = 0.d0
    do j = 1, n
      if (ered .gt. eint(j)) then
        nopen = nopen + 1
      endif
    enddo
    if (nopen .le. 0) then
      ierr = 0
      seps = 0.d0
      nopen_out = 0
      return
    endif
    nopen_out = nopen

    allocate(kop(nopen,nopen), sop(nopen,nopen), qop(nopen,nopen), dk(nopen,nopen), tran(nopen,nopen), phase(nopen))
    kop = (0.d0,0.d0)
    dk = (0.d0,0.d0)
    j = 0
    do jc = 1, n
      if (ered .gt. eint(jc)) then
        j = j + 1
        i = 0
        do ic = 1, n
          if (ered .gt. eint(ic)) then
            i = i + 1
            kop(i,j) = kmat(ic,jc)
            dk(i,j) = dkde(ic,jc)
          endif
        enddo
      endif
    enddo
    call asymptotic_driver_complex(kop, dk, nopen, ederiv, sop, phase, tran, qop, ierr)
    if (ierr .ne. 0) then
      deallocate(kop, sop, qop, dk, tran, phase)
      return
    endif

    ! S-matrix unitarity check: sum of squared moduli of each row should be 1.0
    ! for purely elastic/inelastic scattering. Deviation measures numerical error.
    sumv = 0.d0
    do i = 1, nopen
      sumv = sumv + abs(1.d0 - sum(abs(sop(i,:))**2))
    enddo
    seps = sumv / dble(nopen)

    if (seps .gt. 1.d-2) then
        write(0,*) "WARNING: Large S-matrix unitarity deviation: ", seps, " at E=", ered
    endif

    phase_full(1:nopen) = phase

    if (ederiv) then
      allocate(qwork(nopen,nopen))
      qwork = qop
      ! The collision lifetime matrix Q is defined as Q = i*hbar*S*(dS^dagger/dE).
      ! The eigenvalues of Q are the collision lifetimes [Smith, Phys. Rev. 118, 349 (1960)].
      call hermitian_eig_complex(qwork, nopen, qeval(1:nopen), ierr)
      deallocate(qwork)
      if (ierr .ne. 0) then
        deallocate(kop, sop, qop, dk, tran, phase)
        return
      endif
      do i = 1, nopen
        if (qeval(i) .lt. 0.d0) qeval(i) = -qeval(i)
      enddo
    endif

    j = 0
    do jc = 1, n
      if (ered .gt. eint(jc)) then
        j = j + 1
        i = 0
        do ic = 1, n
          if (ered .gt. eint(ic)) then
            i = i + 1
            cphase = exp(dcmplx(0.d0, -(llev(ic)+llev(jc))*0.5d0*pi))
            smat(ic,jc) = sop(i,j) * cphase
            qphase = exp(dcmplx(0.d0, -(llev(ic)-llev(jc))*0.5d0*pi))
            qmat(ic,jc) = qop(i,j) * qphase
          endif
        enddo
      endif
    enddo
    deallocate(kop, sop, qop, dk, tran, phase)
  end subroutine smatrx_complex

end module mod_smatrx_complex
