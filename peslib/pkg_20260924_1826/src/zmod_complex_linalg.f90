module mod_complex_linalg
  implicit none
  private
  public :: hermitian_eig_complex, solve_linear_complex

  interface
     subroutine zheev(jobz,uplo,n,a,lda,w,work,lwork,rwork,info)
       character(len=1), intent(in) :: jobz, uplo
       integer, intent(in) :: n, lda, lwork
       complex*16, intent(inout) :: a(lda,*)
       double precision, intent(out) :: w(*)
       complex*16, intent(inout) :: work(*)
       double precision, intent(inout) :: rwork(*)
       integer, intent(out) :: info
     end subroutine zheev

     subroutine zgesv(n,nrhs,a,lda,ipiv,b,ldb,info)
       integer, intent(in) :: n, nrhs, lda, ldb
       complex*16, intent(inout) :: a(lda,*)
       integer, intent(out) :: ipiv(*)
       complex*16, intent(inout) :: b(ldb,*)
       integer, intent(out) :: info
     end subroutine zgesv
  end interface

contains

  subroutine hermitian_eig_complex(a, n, w, ierr)
    integer, intent(in) :: n
    complex*16, intent(inout) :: a(n,n)
    double precision, intent(out) :: w(n)
    integer, intent(out) :: ierr

    integer :: lwork
    complex*16, allocatable :: work(:)
    double precision, allocatable :: rwork(:)

    lwork = max(1, 4*n)
    allocate(work(lwork), rwork(max(1, 3*n-2)))
    call zheev('V', 'L', n, a, n, w, work, lwork, rwork, ierr)
    deallocate(work, rwork)
  end subroutine hermitian_eig_complex

  subroutine solve_linear_complex(a, b, n, nrhs, ierr)
    integer, intent(in) :: n, nrhs
    complex*16, intent(inout) :: a(n,n)
    complex*16, intent(inout) :: b(n,nrhs)
    integer, intent(out) :: ierr

    integer, allocatable :: ipiv(:)

    allocate(ipiv(n))
    call zgesv(n, nrhs, a, n, ipiv, b, n, ierr)
    deallocate(ipiv)
  end subroutine solve_linear_complex

end module mod_complex_linalg
