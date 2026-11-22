module bkmp2_c
  use iso_c_binding
  use mod_mkmp2, only: bkmp2_impl => bkmp2
  implicit none
contains
  subroutine pot3_c(r3, v) bind(c, name="pot3_c")
    real(c_double), intent(in)  :: r3(3)
    real(c_double), intent(out) :: v
    real(8) :: dv(3)
    integer :: id
    id = 0
    call bkmp2_impl(r3, v, dv, id)
  end subroutine pot3_c
  subroutine bkmp2_f77(r3, v, dv, id) bind(c, name="bkmp2_")
    real(c_double), intent(in)  :: r3(3)
    real(c_double), intent(out) :: v
    real(8) :: dv(3)
    integer(c_int), intent(in) :: id
    integer :: idd
    idd = id
    call bkmp2_impl(r3, v, dv, idd)
  end subroutine bkmp2_f77
end module bkmp2_c
