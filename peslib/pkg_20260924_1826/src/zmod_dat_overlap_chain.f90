module mod_dat_overlap_chain
  implicit none
  private
  public :: combine_dat_overlaps, transform_logder_complex

contains

  subroutine combine_dat_overlaps(overlaps, n, nsteps, u)
    integer, intent(in) :: n, nsteps
    complex*16, intent(in) :: overlaps(n,n,nsteps)
    complex*16, intent(out) :: u(n,n)
    complex*16 :: tmp(n,n)
    integer :: i, istep

    u = (0.d0, 0.d0)
    do i = 1, n
      u(i,i) = (1.d0, 0.d0)
    enddo
    do istep = 1, nsteps
      call zgemm('N','N',n,n,n,(1.d0,0.d0),u,n,overlaps(:,:,istep),n,(0.d0,0.d0),tmp,n)
      u = tmp
    enddo
  end subroutine combine_dat_overlaps

  subroutine transform_logder_complex(yin, u, n, yout)
    integer, intent(in) :: n
    complex*16, intent(in) :: yin(n,n), u(n,n)
    complex*16, intent(out) :: yout(n,n)
    complex*16 :: tmp(n,n)

    call zgemm('N','N',n,n,n,(1.d0,0.d0),yin,n,u,n,(0.d0,0.d0),tmp,n)
    call zgemm('C','N',n,n,n,(1.d0,0.d0),u,n,tmp,n,(0.d0,0.d0),yout,n)
  end subroutine transform_logder_complex

end module mod_dat_overlap_chain
