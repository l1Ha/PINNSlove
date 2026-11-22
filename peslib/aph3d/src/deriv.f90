!*==DERIV.f90  processed by SPAG 6.55Dc at 15:04 on 21 Aug 2007
SUBROUTINE deriv(t1,t2,t3,c1,c2,c3,pot,vtheta,vchi)                     !      0
USE Numeric_Kinds_Module
USE FileUnits_Module
IMPLICIT NONE
!*--DERIV6
REAL(KIND=WP_Kind) t1, t2, t3, c1, c2, c3, vtheta, vchi, pot(3,3), delta
!ccc   WRITE(Out_Unit,*)'Thetas=',t1,t2,t3
!ccc   WRITE(Out_Unit,*)'Chis=',c1,c2,c3
!ccc   WRITE(Out_Unit,1)(pot(1,j),j=1,3)
!ccc   WRITE(Out_Unit,1)(pot(2,j),j=1,3)
!ccc   WRITE(Out_Unit,1)(pot(3,j),j=1,3)
!ccc 1 FORMAT(3f10.5)
 
delta = t2 - t1                                                         !      0
IF(ABS((delta-(t3-t2))/delta)>1.E-2)THEN
   WRITE(Out_Unit,*)'Error in theta'                                    !      0
   WRITE(Out_Unit,*)t1, t2, t3, t2 - t1, t3 - t2
   STOP
ENDIF
 
VTheta = 0.5*(pot(3,2)-pot(1,2))/delta                                  !      0
 
delta = c2 - c1
IF(ABS((delta-(c3-c2))/delta)>1.E-2)THEN
   WRITE(Out_Unit,*)'Error in chi'                                      !      0
   WRITE(Out_Unit,*)c1, c2, c3, c2 - c1, c3 - c2
   STOP
ENDIF
VChi = 0.5*(pot(2,3)-pot(2,1))/delta                                    !      0
ENDSUBROUTINE DERIV
