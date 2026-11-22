module mod_asymptotic_driver_complex
  use mod_asymptotic_complex, only: symmetrize_hermitian_complex,      &
                                     s_matrix_from_k_complex,           &
                                     s_matrix_derivative_from_k_complex,&
                                     eigenphase_from_k_complex,         &
                                     smith_lifetime_complex
  implicit none
  private
  public :: asymptotic_driver_complex

contains

  subroutine asymptotic_driver_complex(kmat, dkde, nopen, ederiv, smat, phase, phase_tran, qmat, ierr)
    integer, intent(in) :: nopen
    logical, intent(in) :: ederiv
    complex*16, intent(inout) :: kmat(nopen,nopen)
    complex*16, intent(in) :: dkde(nopen,nopen)
    complex*16, intent(out) :: smat(nopen,nopen), phase_tran(nopen,nopen), qmat(nopen,nopen)
    double precision, intent(out) :: phase(nopen)
    integer, intent(out) :: ierr

    complex*16 :: kcopy(nopen,nopen), keig(nopen,nopen)
    complex*16 :: dkcopy(nopen,nopen), dsde(nopen,nopen)

    kcopy = kmat
    call symmetrize_hermitian_complex(kcopy, nopen)
    call s_matrix_from_k_complex(kcopy, nopen, smat, ierr)
    if (ierr .ne. 0) return

    kmat = kcopy
    keig = kcopy
    call eigenphase_from_k_complex(keig, nopen, phase, phase_tran, ierr)
    if (ierr .ne. 0) return

    qmat = (0.d0, 0.d0)
    if (ederiv) then
      dkcopy = dkde
      call symmetrize_hermitian_complex(dkcopy, nopen)
      call s_matrix_derivative_from_k_complex(kcopy, dkcopy, nopen, dsde, ierr)
      if (ierr .ne. 0) return
      call smith_lifetime_complex(smat, dsde, nopen, 27.211386245988d0*2.4188843265857d-2, qmat)
    endif
  end subroutine asymptotic_driver_complex

end module mod_asymptotic_driver_complex
