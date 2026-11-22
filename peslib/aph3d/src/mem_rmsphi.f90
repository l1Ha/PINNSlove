SUBROUTINE mem_rmsphi(nel, numnp, nq, ntheta, nchi, kase)
USE point_rmsphi_Module
USE point_RSF_Module
USE FileUnits_Module
CHARACTER(LEN=*) kase
LOGICAL little, medium, full
INTEGER ithcall, ithsub, ntheta, nchi, nel, numnp, nq

!-----------------------------------------------------------------------
!  Determine printing options.
!-----------------------------------------------------------------------
DATA ithcall/0/, ithsub/0/
DATA little/.false./, medium/.false./, full/.false./
CALL popt ('mem_rmsphi', little, medium, full, ithcall, ithsub)

IF(little)THEN
   WRITE(Out_Unit,*)'Starting to ',kase,' memory in mem_rmsphi'
   WRITE(Out_Unit,*)'      name    length     start    finish'
ENDIF

IF(kase=='ALLOCATE'.or.kase=='release')THEN
   IF(nel>=1)THEN
      ALLOCATE(id(9,nel))
   ELSE
      WRITE(Out_Unit,*)'Error nel must be greater than or equal to 1'
      STOP 'mem_rmsphi'
   ENDIF

   IF(numnp>=9.AND.nq>=1)THEN
      ALLOCATE(th(9,nel))
      ALLOCATE(ch(9,nel))
      ALLOCATE(phi(numnp,nel))
   ELSE
      WRITE(Out_Unit,*)'Error numnp and nq must be greater than or','equal to 9'
      WRITE(Out_Unit,*)'numnp = ',numnp,'  nq = ',nq
      STOP 'mem_rmsphi'
   ENDIF

   IF(ntheta>=1.AND.nchi>1)THEN
      ALLOCATE(phirms(ntheta,nchi))
   ELSE
      WRITE(Out_Unit,*)'Error ntheta and nchi must be greater than or',' equal to 1'
   ENDIF

   IF(little)THEN
   WRITE(Out_Unit,*)'Memory ',kase,'d in mem_rmsphi'
   ENDIF

ELSE
   WRITE(Out_Unit,*)'Error kase /= ALLOCATE or release'
   STOP 'mem_rmsphi'
ENDIF

RETURN
END
