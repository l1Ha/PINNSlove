module mod_scratchfile_complex
  implicit none
  private
  public :: setrec_complex, putrec_complex, getrec_complex, endrec_complex

contains

  subroutine putrec_complex(a, n, nrg)
    integer, intent(in) :: n, nrg
    complex*16, intent(in) :: a(n,n)
    integer :: irec

    irec = nrg + 1
    write(unit=11, rec=irec) a
  end subroutine putrec_complex

  subroutine setrec_complex(n, nnrg)
    integer, intent(in) :: n, nnrg
    double precision :: diskmb

    diskmb = 1.6d-5 * n * n * (nnrg + 1)
    open(unit=11, status='scratch', form='unformatted', access='direct', recl=16*n*n)
    write(6,61) diskmb
61  format(/1x,'SETREC_COMPLEX:'/1x,70('-')/1x, &
         '*** This run is using ',f8.3,'  Mb of complex scratch disk space'/1x, &
         70('-'))
  end subroutine setrec_complex

  subroutine endrec_complex
    close(unit=11, status='delete')
  end subroutine endrec_complex

  subroutine getrec_complex(a, n, nrg)
    integer, intent(in) :: n, nrg
    complex*16, intent(out) :: a(n,n)
    integer :: irec

    irec = nrg + 1
    read(unit=11, rec=irec) a
  end subroutine getrec_complex

end module mod_scratchfile_complex
