SUBROUTINE mem_quad (nquad, ntheta, nchi, nq)
USE FileUnits_Module
USE point_quad_Module
CHARACTER(LEN=8) kase
LOGICAL little, medium, full
INTEGER ithcall, ithsub, nquad, ntheta, nchi, nq
DATA kase/'ALLOCATE'/

!----------------------------------------------------------------------
!  Determine printing options.
!----------------------------------------------------------------------
DATA ithcall/0/, ithsub/0/
DATA little/.false./, medium/.false./, full/.false./
CALL popt ('mem_quad', little, medium, full, ithcall, ithsub)

!----------------------------------------------------------------------
! Allocate memory for quadrature arrays.
!----------------------------------------------------------------------
IF(little)THEN
   WRITE(Out_unit,*)'Allocating memory in mem_quad'
   WRITE(Out_unit,*)'      name    length     start    finish'
ENDIF


IF(nquad>=3)THEN
   ALLOCATE(xpt(nquad))
   ALLOCATE(wht(nquad))
   ALLOCATE(value(nquad,nquad))
   ALLOCATE(h(9,nquad,nquad))
   ALLOCATE(p(2,9,nquad,nquad))
   ALLOCATE(hcoef(9,nquad,nquad))
ELSE
   WRITE(Out_unit,*)'Error nquad must be greater than or equal to 3'
   STOP 'mem_quad'
ENDIF

IF(ntheta>0)THEN
   ALLOCATE(thetaval(ntheta*nchi))
ELSE
   WRITE(Out_unit,*)'Error ntheta must be greater than or equal to 1'
   STOP 'mem_quad'
ENDIF

IF(nchi>0)THEN
   ALLOCATE(chivals(ntheta*nchi))
ELSE
   WRITE(Out_unit,*)'Error nchi must be greater than or equal to 1 '
   STOP 'mem_quad'
ENDIF

IF(nq>0)THEN
   ALLOCATE(ovrlapmx(nq))
   ALLOCATE(energy(nq))
   ALLOCATE(energy1(nq))
   ALLOCATE(lambda(nq))
ELSE
   WRITE(Out_unit,*)'Error nq must be greater than or equal to 1'
   STOP 'mem_quad'

ENDIF

IF(little)THEN
WRITE(Out_unit,*)'Memory allocated in mem_quad'
ENDIF

RETURN
END
