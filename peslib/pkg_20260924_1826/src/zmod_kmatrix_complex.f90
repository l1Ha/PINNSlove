module mod_kmatrix_complex
  use mod_complex_linalg, only: solve_linear_complex
  implicit none
  private
  public :: kmatrix_from_r_complex, kmatrix_derivative_from_r_complex

contains

  subroutine kmatrix_from_r_complex(rmat, amat, bmat, emat, fmat, n, kmat, ierr)
    integer, intent(in) :: n
    complex*16, intent(in) :: rmat(n,n), amat(n,n), bmat(n,n), emat(n,n), fmat(n,n)
    complex*16, intent(out) :: kmat(n,n)
    integer, intent(out) :: ierr

    complex*16 :: lhs(n,n), rhs(n,n)

    lhs = -(matmul(rmat, fmat) - bmat)
    rhs =  (matmul(rmat, emat) - amat)
    call solve_linear_complex(lhs, rhs, n, n, ierr)
    kmat = rhs
  end subroutine kmatrix_from_r_complex

  subroutine kmatrix_derivative_from_r_complex(rmat, drmat, amat, bmat, emat, fmat, damat, dbmat, demat, dfmat, n, kmat, dkde, ierr)
    integer, intent(in) :: n
    complex*16, intent(in) :: rmat(n,n), drmat(n,n)
    complex*16, intent(in) :: amat(n,n), bmat(n,n), emat(n,n), fmat(n,n)
    complex*16, intent(in) :: damat(n,n), dbmat(n,n), demat(n,n), dfmat(n,n)
    complex*16, intent(in) :: kmat(n,n)
    complex*16, intent(out) :: dkde(n,n)
    integer, intent(out) :: ierr

    complex*16 :: lhs(n,n), rhs(n,n), tmp(n,n)

    lhs = -(matmul(rmat, fmat) - bmat)
    tmp = matmul(drmat, fmat) + matmul(rmat, dfmat) - dbmat
    rhs = matmul(drmat, emat) + matmul(rmat, demat) - damat - matmul(tmp, kmat)
    call solve_linear_complex(lhs, rhs, n, n, ierr)
    dkde = rhs
  end subroutine kmatrix_derivative_from_r_complex

end module mod_kmatrix_complex
