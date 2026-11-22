module mod_exchange_demix
  implicit none
  private
  public :: demix_exchange_s

contains

  subroutine demix_exchange_s(sdir, sloop, seven, sodd)
    complex*16, intent(in) :: sdir, sloop
    complex*16, intent(out) :: seven, sodd
    complex*16 :: fdir, flo

    fdir = sdir - (1.d0, 0.d0)
    flo = sloop
    seven = (1.d0, 0.d0) + (fdir + flo) / sqrt(2.d0)
    sodd  = (1.d0, 0.d0) + (fdir - flo) / sqrt(2.d0)
  end subroutine demix_exchange_s

end module mod_exchange_demix
