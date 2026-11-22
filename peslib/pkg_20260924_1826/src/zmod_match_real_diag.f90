module mod_match_real_diag
  private
  public :: build_k_matrix_real

contains

  subroutine build_k_matrix_real(nrg, kmat, cvr, cvi, cro, eint, ilev, jlev, &
                                 klev, nlev, ered)
    use mod_common
    use Math_Functions, only: bessel
    use mod_scratchfile, only: getrec
    use mod_frames, only: frames
    implicit double precision (a-h,o-z)

    integer, intent(in) :: nrg
    double precision, intent(out) :: kmat(nchan,nchan)
    double precision, intent(in) :: cvr(nvi,n), cvi(nvi,n), cro(3,n)
    double precision, intent(in) :: eint(n), ered(*)
    integer, intent(in) :: ilev(n), jlev(n), klev(n), nlev(n)

    double precision :: y(nchan,nchan), z(nchan,nchan), cmat(nchan,nchan)
    double precision :: llev(nchan)

    call getrec(y, nchan, nrg)
    call frames(y, cmat, cro, jlev, klev, llev, +1)
    call bessel(kmat, y, z, cvr, cvi, eint, ilev, jlev, klev, llev, nlev, &
                ered(nrg))
  end subroutine build_k_matrix_real

end module mod_match_real_diag
