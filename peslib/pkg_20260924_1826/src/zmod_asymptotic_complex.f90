module mod_asymptotic_complex
  use mod_complex_linalg, only: hermitian_eig_complex
  implicit none
  private
  public :: symmetrize_hermitian_complex
  public :: s_matrix_from_k_complex, s_matrix_derivative_from_k_complex
  public :: eigenphase_from_k_complex
  public :: smith_lifetime_complex

  interface
    subroutine zgemm(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
      character(len=1), intent(in) :: transa, transb
      integer, intent(in) :: m, n, k, lda, ldb, ldc
      complex*16, intent(in) :: alpha, beta
      complex*16, intent(in) :: a(lda,*), b(ldb,*)
      complex*16, intent(inout) :: c(ldc,*)
    end subroutine zgemm
  end interface

contains

  subroutine symmetrize_hermitian_complex(a, n)
    integer, intent(in) :: n
    complex*16, intent(inout) :: a(n,n)
    integer :: i, j
    complex*16 :: tmp

    do j = 1, n
       do i = 1, j-1
          tmp = 0.5d0 * (a(i,j) + dconjg(a(j,i)))
          a(i,j) = tmp
          a(j,i) = dconjg(tmp)
       enddo
       a(j,j) = dcmplx(dble(a(j,j)), 0.d0)
    enddo
  end subroutine symmetrize_hermitian_complex

  subroutine s_matrix_from_k_complex(kmat, n, smat, ierr)
    ! Compute S = (I - iK)^{-1}(I + iK) via eigendecomposition of the
    ! Hermitian K-matrix.  This is numerically stable even when ||K|| >> 1,
    ! whereas the (I + K^2) direct-solve formulation is ill-conditioned for
    ! large K-matrix eigenvalues.
    integer, intent(in) :: n
    complex*16, intent(in) :: kmat(n,n)
    complex*16, intent(out) :: smat(n,n)
    integer, intent(out) :: ierr

    complex*16, allocatable :: kcopy(:,:), tmp(:,:)
    double precision, allocatable :: evals(:)
    integer :: i, j
    complex*16 :: phase

    allocate(kcopy(n,n), evals(n), tmp(n,n))

    kcopy = kmat
    call symmetrize_hermitian_complex(kcopy, n)
    ! kcopy returns as the eigenvector matrix U after the call
    call hermitian_eig_complex(kcopy, n, evals, ierr)
    if (ierr .ne. 0) then
      deallocate(kcopy, evals, tmp)
      return
    endif

    ! S_eig(i) = (1 + i*k_i) / (1 - i*k_i)  for each eigenvalue k_i
    ! Build tmp = U * diag(S_eig) * U^dagger
    do j = 1, n
      phase = dcmplx(1.d0, evals(j)) / dcmplx(1.d0, -evals(j))
      do i = 1, n
        tmp(i,j) = kcopy(i,j) * phase
      enddo
    enddo
    ! smat = tmp * U^dagger = tmp * conjg(transpose(kcopy_orig))
    call zgemm('N','C',n,n,n,(1.d0,0.d0),tmp,n,kcopy,n,(0.d0,0.d0),smat,n)

    deallocate(kcopy, evals, tmp)
  end subroutine s_matrix_from_k_complex

  subroutine s_matrix_derivative_from_k_complex(kmat, dkde, n, dsde, ierr)
    ! Compute dS/dE = 2i (I-iK)^{-1} dK/dE (I-iK)^{-1} using the
    ! eigendecomposition of the Hermitian K-matrix for numerical stability.
    integer, intent(in) :: n
    complex*16, intent(in) :: kmat(n,n), dkde(n,n)
    complex*16, intent(out) :: dsde(n,n)
    integer, intent(out) :: ierr

    complex*16, allocatable :: kcopy(:,:), dkrot(:,:), tmp(:,:)
    double precision, allocatable :: evals(:)
    integer :: i, j

    allocate(kcopy(n,n), evals(n), dkrot(n,n), tmp(n,n))

    kcopy = kmat
    call symmetrize_hermitian_complex(kcopy, n)
    ! After call, kcopy holds the eigenvector matrix U; evals holds k_i
    call hermitian_eig_complex(kcopy, n, evals, ierr)
    if (ierr .ne. 0) then
      deallocate(kcopy, evals, dkrot, tmp)
      return
    endif

    ! dkrot = U^dagger * dK/dE * U  (rotate derivative to eigenbasis)
    call zgemm('C','N',n,n,n,(1.d0,0.d0),kcopy,n,dkde,n,(0.d0,0.d0),tmp,n)
    call zgemm('N','N',n,n,n,(1.d0,0.d0),tmp,n,kcopy,n,(0.d0,0.d0),dkrot,n)

    ! Apply diagonal factors: dkrot(i,j) *= 2i / ((1 - i*k_i)(1 - i*k_j))
    do j = 1, n
      do i = 1, n
        dkrot(i,j) = (0.d0,2.d0) * dkrot(i,j) &
                     / (dcmplx(1.d0,-evals(i)) * dcmplx(1.d0,-evals(j)))
      enddo
    enddo

    ! Transform back: dsde = U * dkrot * U^dagger
    call zgemm('N','N',n,n,n,(1.d0,0.d0),kcopy,n,dkrot,n,(0.d0,0.d0),tmp,n)
    call zgemm('N','C',n,n,n,(1.d0,0.d0),tmp,n,kcopy,n,(0.d0,0.d0),dsde,n)

    deallocate(kcopy, evals, dkrot, tmp)
  end subroutine s_matrix_derivative_from_k_complex

  subroutine eigenphase_from_k_complex(kmat, n, phase, tran, ierr)
    integer, intent(in) :: n
    complex*16, intent(inout) :: kmat(n,n)
    double precision, intent(out) :: phase(n)
    complex*16, intent(out) :: tran(n,n)
    integer, intent(out) :: ierr

    double precision :: evals(n)
    integer :: i

    call symmetrize_hermitian_complex(kmat, n)
    call hermitian_eig_complex(kmat, n, evals, ierr)
    if (ierr .ne. 0) return
    tran = kmat
    do i = 1, n
      phase(i) = atan(evals(i))
    enddo
  end subroutine eigenphase_from_k_complex

  subroutine smith_lifetime_complex(smat, dsde, n, hbar_ev, qmat)
    integer, intent(in) :: n
    complex*16, intent(in) :: smat(n,n), dsde(n,n)
    double precision, intent(in) :: hbar_ev
    complex*16, intent(out) :: qmat(n,n)
    complex*16 :: tmp(n,n)

    call zgemm('N','C',n,n,n,(1.d0,0.d0),dsde,n,smat,n,(0.d0,0.d0),tmp,n)
    qmat = -(0.d0,1.d0) * hbar_ev * tmp
    call symmetrize_hermitian_complex(qmat, n)
  end subroutine smith_lifetime_complex

end module mod_asymptotic_complex
