module mod_asymptotic_real
  use mod_symevp, only: symevp
  use mod_complex_linalg, only: hermitian_eig_complex
  implicit none
  private
  public :: asymptotic_real_diagnostics

  interface
     subroutine dgesv(n,nrhs,a,lda,ipiv,b,ldb,info)
       integer, intent(in) :: n, nrhs, lda, ldb
       double precision, intent(inout) :: a(lda,*), b(ldb,*)
       integer, intent(out) :: ipiv(*)
       integer, intent(out) :: info
     end subroutine dgesv
  end interface

contains

  subroutine solve_linear_real(a, b, n, nrhs, ierr)
    integer, intent(in) :: n, nrhs
    double precision, intent(inout) :: a(n,n), b(n,nrhs)
    integer, intent(out) :: ierr

    integer, allocatable :: ipiv(:)

    allocate(ipiv(n))
    call dgesv(n, nrhs, a, n, ipiv, b, n, ierr)
    deallocate(ipiv)
  end subroutine solve_linear_real

  subroutine asymptotic_real_diagnostics(kmat, dkde, eint, n, ered, ederiv, &
                                         phase_full, qdiag_full, qeval_full, &
                                         nopen_out, ierr)
    integer, intent(in) :: n
    double precision, intent(in) :: kmat(n,n), dkde(n,n), eint(n), ered
    logical, intent(in) :: ederiv
    double precision, intent(out) :: phase_full(n), qdiag_full(n), qeval_full(n)
    integer, intent(out) :: nopen_out, ierr

    double precision, allocatable :: kop(:,:), dkop(:,:), kdiag(:,:), evals(:)
    double precision, allocatable :: unitmat(:,:), sreal(:,:), simag(:,:)
    double precision, allocatable :: se_real(:,:), se_imag(:,:), rhs(:,:)
    double precision, allocatable :: qreal(:,:), qimag(:,:)
    complex*16, allocatable :: qmat(:,:)
    integer :: i, j, ic, jc, nopen
    double precision, parameter :: hbar_ev_fs = 27.211386245988d0 * &
                                               2.4188843265857d-2

    phase_full = 0.d0
    qdiag_full = 0.d0
    qeval_full = 0.d0
    ierr = 0

    nopen = 0
    do i = 1, n
      if (ered .gt. eint(i)) nopen = nopen + 1
    enddo
    nopen_out = nopen
    if (nopen .le. 0) return

    allocate(kop(nopen,nopen), dkop(nopen,nopen), kdiag(nopen,nopen), evals(nopen))

    jc = 0
    do j = 1, n
      if (ered .gt. eint(j)) then
        jc = jc + 1
        ic = 0
        do i = 1, n
          if (ered .gt. eint(i)) then
            ic = ic + 1
            kop(ic,jc) = 0.5d0 * (kmat(i,j) + kmat(j,i))
            dkop(ic,jc) = 0.5d0 * (dkde(i,j) + dkde(j,i))
          endif
        enddo
      endif
    enddo

    kdiag = kop
    call symevp(kdiag, nopen, nopen, evals, ierr)
    if (ierr .ne. 0) then
      deallocate(kop, dkop, kdiag, evals)
      return
    endif
    do i = 1, nopen
      phase_full(i) = atan(evals(i))
    enddo

    if (.not. ederiv) then
      deallocate(kop, dkop, kdiag, evals)
      return
    endif

    allocate(unitmat(nopen,nopen), sreal(nopen,nopen), simag(nopen,nopen))
    allocate(se_real(nopen,nopen), se_imag(nopen,nopen), rhs(nopen,nopen))
    allocate(qreal(nopen,nopen), qimag(nopen,nopen), qmat(nopen,nopen))

    unitmat = 0.d0
    do i = 1, nopen
      unitmat(i,i) = 1.d0
    enddo

    kdiag = unitmat + matmul(kop, kop)
    rhs = 2.d0 * kop
    call solve_linear_real(kdiag, rhs, nopen, nopen, ierr)
    if (ierr .ne. 0) then
      deallocate(kop, dkop, kdiag, evals, unitmat, sreal, simag)
      deallocate(se_real, se_imag, rhs, qreal, qimag, qmat)
      return
    endif
    simag = rhs
    sreal = unitmat - matmul(simag, kop)

    kdiag = unitmat + matmul(kop, kop)
    rhs = -matmul(matmul(dkop, kop) + matmul(kop, dkop), simag) + 2.d0 * dkop
    call solve_linear_real(kdiag, rhs, nopen, nopen, ierr)
    if (ierr .ne. 0) then
      deallocate(kop, dkop, kdiag, evals, unitmat, sreal, simag)
      deallocate(se_real, se_imag, rhs, qreal, qimag, qmat)
      return
    endif
    se_imag = rhs
    se_real = -matmul(se_imag, kop) - matmul(simag, dkop)

    qreal = hbar_ev_fs * (matmul(se_imag, transpose(sreal)) - &
             matmul(se_real, transpose(simag)))
    qimag = -hbar_ev_fs * (matmul(se_real, transpose(sreal)) + &
             matmul(se_imag, transpose(simag)))

    do i = 1, nopen
      qdiag_full(i) = qreal(i,i)
      do j = 1, nopen
        qmat(i,j) = dcmplx(qreal(i,j), qimag(i,j))
      enddo
    enddo
    do j = 1, nopen
      do i = 1, j
        qmat(i,j) = 0.5d0 * (qmat(i,j) + dconjg(qmat(j,i)))
        qmat(j,i) = dconjg(qmat(i,j))
      enddo
    enddo

    call hermitian_eig_complex(qmat, nopen, qeval_full(1:nopen), ierr)

    deallocate(kop, dkop, kdiag, evals, unitmat, sreal, simag)
    deallocate(se_real, se_imag, rhs, qreal, qimag, qmat)
  end subroutine asymptotic_real_diagnostics

end module mod_asymptotic_real
