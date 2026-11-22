subroutine testpot(jval,mu)
USE Numeric_Kinds_Module
USE FileUnits_Module
IMPLICIT NONE

integer :: i, j
integer :: jval

REAL(KIND=WP_Kind) :: mu, hstep

REAL(KIND=WP_Kind) :: r(23,2),vrkr(23),v(23,2), vpp, vp


open(unit=Table4_dat_Unit,file="table4.dat")
do i=1,23
    read(Table4_dat_Unit,*) r(i,1), r(i,2), vrkr(i)
enddo

r=r/0.5291772108d0

do i=1,23
    CALL PotDiatomic(r(i,1), hstep, v(i,1),vp, vpp, jval, mu)
    CALL PotDiatomic(r(i,2), hstep, v(i,2),vp, vpp, jval, mu)
enddo

v=v*2.194746313705d5

do i=1,23
    v(i,1)=v(i,1)-v(23,1)
    v(i,2)=v(i,2)-v(23,2)
enddo

open(unit=TestPot_Unit,file="tab4comp.dat")
write(TestPot_Unit,*)mu
do i=1,23
    write(TestPot_Unit,*) r(i,1), vrkr(i), v(i,1)
enddo
do j=1,23
    i=24-j
    write(TestPot_Unit,*) r(i,2), vrkr(i), v(i,2)
enddo

end subroutine testpot

